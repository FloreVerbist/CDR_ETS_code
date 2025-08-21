function ADMM_price_based_m!(Case_name::String, results::Dict,ADMM::Dict, model_e::Model, model_r::Model, model_ccb::Model, dict::Dict, ETS_data::Dict, LINEAR::Bool) #, conditional_case::Bool)
    convergence = 0
    iterations = ProgressBar(1:ADMM_max_it-1)
    Cum_remov = model_r.ext[:parameters][:Cum_remov] 
    Cum_abat = model_e.ext[:parameters][:Cum_abat] 
    Beta = model_r.ext[:parameters][:Beta] 
    Gamma = model_r.ext[:parameters][:Gamma] 


    Emitter!(model_e)
    

    #AC(abat, t) = @. (D_a ./ (1 + Cum_abat[t-1]).^(M_a) .*(abat./((U_a/(1+exp(K_a*(t+1)))) .- abat)).^C_a)

    RC(remov, t, Cum_remov) = f(remov, t, Cum_remov) # NOTE: this is different compared to all other cases
    AC(abat, t, Cum_abat) = f2(abat, t, Cum_abat) # NOTE: this is different compared to all other cases

    for iter in iterations
        if convergence == 0
            # Calculate penalty terms ADMM and update price to most recent value
            if iter > 1 # defining previous iteration values 

  
                model_e.ext[:parameters][:EAC_price] = Dict(t => results["EAC_price"][iter,t+1] for t in model_e.ext[:sets][:T])
                model_e.ext[:parameters][:ec_bar] = Dict(t => results["ec"][iter-1,t+1] for t in model_e.ext[:sets][:T])



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
                
        
                D_r = model_r.ext[:parameters][:D_r]  
                C_r = model_r.ext[:parameters][:C_r]  
                M_r = model_r.ext[:parameters][:M_r] 
                K_r = model_r.ext[:parameters][:K_r] 
                U_r = model_r.ext[:parameters][:U_r] 
    



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


                abat = collect(0:st:ab_max)
                EAC_price_nom = maximum([0,model_e.ext[:parameters][:EAC_price][t]]) # FV: maximum maybe not correct here?
                LAMDA_ER = @.(EAC_price_nom)
                _, j = findmin(abs.(AC(abat,t, Cum_abat).-LAMDA_ER)) #FV_add: t+1

                Abat_t = maximum([0, abat[j]]).*EM_REF

                model_e.ext[:parameters][:e][t] = EM_REF - maximum([0,Abat_t])

                Cum_abat[t] = model_e.ext[:parameters][:Cum_abat][t] = model_e.ext[:parameters][:Cum_abat][t-1] + Abat_t
                
                results["Cum_abat"][iter, t+1] = Cum_abat[t]

                # model_e.ext[:parameters][:e][t] = maximum([(EM_REF - (LAMDA_ER/(BETA)).^(1/GAMMA)), 0])

                #----------------------------

            end

            results["e"][iter,:] = [model_e.ext[:parameters][:e][t] for t in model_e.ext[:sets][:T]]
 
            Update_emitter!(model_e)
            variables_extract(model_e)

            results["ec"][iter,:] = [value.(model_e.ext[:variables][:ec][t]) for t in model_e.ext[:sets][:T]]


              
            # # Imbalances
            ADMM["Imbalances"]["EAC"][iter,:] = [model_e.ext[:parameters][:CAP][t] for t in  model_e.ext[:sets][:T]]  - results["ec"][iter,:] 


            # # Primal residuals
            ADMM["Residuals"]["Primal_EAC"][iter] = sqrt(sum(ADMM["Imbalances"]["EAC"][iter,:].^2))

            if iter > 1 
                ADMM["Residuals"]["Dual_EC"][iter] = sqrt(sum(( (results["ec"][iter,:] - results["ec"][iter-1,:])).^2))
            else
            end
           
            # Price updates 
            results["EAC_price"][iter+1,:] = results["EAC_price"][iter,:] - model_e.ext[:parameters][:ρ] .*ADMM["Imbalances"]["EAC"][iter,:]

            # Progress bar
            set_description(iterations, string(@sprintf("Primal EAC %.3f -- Dual EC %.3f",  ADMM["Residuals"]["Primal_EAC"][iter],ADMM["Residuals"]["Dual_EC"][iter])))

            # Check convergence: primal and dual satisfy tolerance OR primal remains high and dual near zero. 
            # In case the latter happens, the solution needs to be inspected to ensure this can be interpreted as an equilibrium.
            if (ADMM["Residuals"]["Primal_EAC"][iter] <= model_e.ext[:parameters][:epsilon] && ADMM["Residuals"]["Dual_EC"][iter] <= model_e.ext[:parameters][:epsilon])
                convergence = 1
            end
            ADMM["n_iter"] = copy(iter)
        end
    end
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
    CDR_price = model_r.ext[:parameters][:CDR_price] 

    # parameters_extract(model_r)
    # Total_costs_rem = sum(1/(1+DRc)^t*Area(fr, t, Cum_remov_extract, 0, (r[t]./EM_REF+0.00001), 50000).*EM_REF for t in model_r.ext[:sets][:T])[1]
    # Total_costs_rem_NL = sum(1/(1+DRc)^t*Area(fr, t, Cum_NoLearning_remov, 0, (r[t]./EM_REF+0.00001), 50000).*EM_REF for t in model_r.ext[:sets][:T])[1]
    # parameters_extract(model_e)
    # Total_costs_abat = sum(1/(1+DRc)^t*Area(f2, t, Cum_abat_extract, 0, ((EM_REF-e[t])./EM_REF+0.00001), 50000).*EM_REF for t in model_r.ext[:sets][:T])[1]
    # Total_costs_abat_NL = sum(1/(1+DRc)^t*Area(fa, t, Cum_NoLearning_abat, 0, ((EM_REF-e[t])./EM_REF+0.00001), 50000).*EM_REF for t in model_r.ext[:sets][:T])[1]
    # EAC_price = model_e.ext[:parameters][:EAC_price]
    # CDR_price = model_r.ext[:parameters][:CDR_price] # FV: maybe this one should be changed as well. 


    Total_costs = Total_costs_rem + Total_costs_abat
    
   
    Total_budget = sum((results["ec"][ADMM["n_iter"],:].*[values(EAC_price[t]) for t in 0:1:T_OPT])) 
    Total_net_emissions = sum(results["e"][ADMM["n_iter"],:] - results["r"][ADMM["n_iter"],:])
    Total_emissions  = sum(results["e"][ADMM["n_iter"],:]) 
    Total_removals = sum(results["r"][ADMM["n_iter"],:]) 
    Total_abated =  (T_OPT+1)*EM_REF - Total_emissions
    ETS_max = maximum([EAC_price[t] for t in 0:1:T_OPT])
    RC_max = maximum([results["CDR_price"][ADMM["n_iter"],t+1] for t in 0:1:T_OPT])
    ETS_avg = sum([EAC_price[t].*results["ec"][ADMM["n_iter"],t+1] for t in 0:1:T_OPT])./sum([results["ec"][ADMM["n_iter"],t+1] for t in 0:1:T_OPT])
    credit_costs_emitter_avg = sum([EAC_price[t].*results["ec"][ADMM["n_iter"],t+1] for t in 0:1:T_OPT])./sum([results["e"][ADMM["n_iter"],t+1] for t in 0:1:T_OPT])


    push!(dict[(Case_name,"Total_net_emissions")], Total_net_emissions/1000) #Gton
    push!(dict[(Case_name,"Total_emissions")], Total_emissions/1000) #Gton
    push!(dict[(Case_name,"Fraction_removals")], Total_net_emissions/Total_emissions)
    push!(dict[(Case_name,"Total_costs_per_net")], Total_costs/Total_net_emissions)
    push!(dict[(Case_name,"Total_costs_per_abat")], Total_costs/Total_abated)
    push!(dict[(Case_name,"Total_budget_per_net")], Total_budget/Total_net_emissions)
    push!(dict[(Case_name,"Learning_gains_per_rem")], 0.0)
    push!(dict[(Case_name,"Learning_gains_per_em")], (Total_costs_abat_NL- Total_costs_abat)/Total_abated)
    push!(dict[(Case_name,"ETS_max")], ETS_max)
    push!(dict[(Case_name,"ETS_avg")], ETS_avg)
    push!(dict[(Case_name,"RC_max")], RC_max)
    push!(dict[(Case_name,"ETS_avg")], ETS_avg)
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

   


    # Define variables
    #e = model.ext[:variables][:e] = @variable(model, e[T], lower_bound=0, base_name="Emissions bought in year t by industrial fringe")
    ec = model.ext[:variables][:ec] = @variable(model, ec[T], lower_bound=0, base_name="Emission allowances bought in year t by industrial fringe")

    # Constraints 
    Banking_e =  model.ext[:constraints][:Banking_e] = @constraint(model, [t=T],
    (sum(ec[i] for  i in 0:t)  >= sum(e[i]   for i in 0:t)))
    

    # Objective 
    
    model.ext[:objective] = @objective(model, Min, 
    # MACC
    sum(DiscountFactor[t]*(CAP[t].*EAC_price[t]) for t in T) +  sum(DiscountFactor[t]*((ec[t]-CAP[t]).*EAC_price[t]) for t in T) + #sum(DiscountFactor[t]*BETA*(EM_REF - e[t]).^(2)/2 for t in T) +
    # Penalty term ETS market
    ( ρ./2) * sum( DiscountFactor[t]*((ec[t] - (ec_bar[t] + 1/N*(CAP[t] -  ec[t])))^2) for t in T) .* a  )
    
    
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
    sum(DiscountFactor[t]*(CAP[t].*EAC_price[t]) for t in T) +  sum(DiscountFactor[t]*((ec[t]-CAP[t]).*EAC_price[t]) for t in T) + # sum(DiscountFactor[t]*BETA*(EM_REF - e[t]).^(2)/2 for t in T) +
    # Penalty term ETS market
    ( ρ./2) * sum( DiscountFactor[t]*((ec[t] - (ec_bar[t] + 1/N*(CAP[t]  -  ec[t])))^2) for t in T) .* a  )
    
    
    optimize!(model);

    return model 
end


function plot_model!(Case_name::String, C_nr, model_e, model_r, model_ccb, ADMM, results, Dict_results::Dict, plotting::Bool)
    
    if !plotting
        skip
    else
        ymax = 650.0
        parameters_extract(model_e)



        # Step 6: 
        EM = results["e"][ADMM["n_iter"],:]
        hatch_vect = [x for p in zip(EM, zeros(length(EM))) for x in p].*[0; repeat([1, 1, 0, 0], convert(Int64, T_OPT/2)); 1]
        Doubled_time_vector =  [x for x in Timeseries_data[1:1:(T_OPT+1),"Year"] for _ in 1:2] + [repeat([0, 0, 0, -1], convert(Int64, T_OPT/2)); 0; 0]


        Peac = results["EAC_price"][ADMM["n_iter"],:]


        variables_extract(model_e)

        # parameters_extract(model)
        # variables_extract(model)


        EAC_surplus = [sum( ec[i]  - (EM[i+1])  for i in 0:t) for t in model_e.ext[:sets][:T]]

        Color = ["snow3" "snow2"  "slategray3"]
        Alpha = [1 0.5 0.5]

        pp = stackedarea(Timeseries_data[1:1:(T_OPT+1),"Year"], [EM[t+1] for t in 0:1:(T_OPT)],  ylims = (-200,3000), labels = "Emissions", ylabel = L"Emissions [MtCO$_2$/yr]", xlabel = "Years", legend = (0.02, 0.9), foreground_color_legend = nothing, linecolor = Color, color = Color, alpha = Alpha)     # legend = (0.32, -0.12)
        plot!(Doubled_time_vector, [hatch_vect], ylims = (-200, 3000), labels = "", title = "",  linecolor = :white, alpha = 0.4, lw=1)  
        plot!(twiny(pp), Timeseries_data[1:1:(T_OPT+1),"Year"], [Timeseries_data[1:1:(T_OPT+1), "ETS_cap [EUA/y]"]./1], ylims = (-200,3000), label = ["ETS cap" ], lw =2, xformatter=_->"", yformatter=_->"", xlabel = "", ylabel="", legend = (0.02,0.98), foreground_color_legend = nothing, color = ["black"  "navy" "black" "brown"], ls = [:solid :dash :dash :dash])
        # plot!(twiny(pp), Timeseries_data[1:1:(T_OPT+1),"Year"], [Timeseries_data[1:1:(T_OPT+1), "ETS_cap [EUA/y]"]./1 EAC_surplus], ylims = (-200,3000), label = ["EAC cap" "EAC surplus emitter" ], xformatter=_->"", yformatter=_->"", xlabel = "", ylabel="", legend = (0.02,0.98), color = ["black"  "navy" "black" "brown"], ls = [:solid :dash :dash :dash])
        plot!(twinx(pp), Timeseries_data[1:1:(T_OPT+1),"Year"], [[Peac[t]   for t in 1:1:(T_OPT+1)] ], ylims = (0, ymax), yticks = 50:50:ymax, label = ["ETS price (MAC)"], color = ["grey36" "grey83"], ylabel = [L"Certificate price [€/tCO$_2$]"], legend = (0.58,0.98), foreground_color_legend = nothing, lw = 2)
        pprint = plot(pp, size  = (500, 400), xlim = [2024,(2024+T_OPT)],  titlefont = font(11),  title = Plot_title[C_nr], plot_titlefontsize = 12, topmargin = 2mm)
        annotate!(2032, 2570,text("/////", :white, :right, 6), linewidth =2)

        ppp = groupedbar(Timeseries_data[1:1:(T_OPT+1),"Year"], [EM zeros(length(EM))], bar_position =:stack, xticks = (2025:5:(2024+T_OPT)), ylims = (-600, 2500), yticks =-500:250:2500, labels = ["Emissions"  "Removals"], foreground_color_legend = nothing, title = "", ylabel = L"Emissions or Removals [MtCO$_2$/yr]", xlabel = "Years", legend = (0.02, 0.85), lw = 0.75, color = Color, linecolor = [:gray31 :gray31], alpha = Alpha, right_margin = 5mm)
        plot!(twiny(ppp), Timeseries_data[1:1:(T_OPT+1),"Year"], [Timeseries_data[1:1:(T_OPT+1), "ETS_cap [EUA/y]"]./1], lw=[1 1], xticks = (2025:5:(2024+T_OPT)), ylims = (-600, 2500), label = ["ETS cap"],foreground_color_legend = nothing, xformatter=_->"", yformatter=_->"", xlabel = "", ylabel="", legend = (0.02,0.98), color = ["black"  "black" "black" "brown"], ls = [:solid :dash :dash :dash], margin = 0mm)
        plot!(twiny(ppp), Timeseries_data[1:1:(T_OPT+1),"Year"], [EM], lw=0, xticks = (2025:5:(2024+T_OPT)), ylims = (-600, 2500), label = ["Net emissions"],foreground_color_legend = nothing, xformatter=_->"", yformatter=_->"", xlabel = "", ylabel="", legend = (0.02,0.91), color = ["black"  "black" "black" "brown"], ls = [:dash :dash :dash :dash], margin = 0mm, markershape = :xcross, markersize = 2)
        plot!(twinx(ppp), Timeseries_data[1:1:(T_OPT+1),"Year"], [[Peac[t]   for t in 1:1:(T_OPT+1)]], ylims = (0, ymax), yticks = 0:50:ymax, label = ["ETS price (MAC)" "RC price (MRC)"], foreground_color_legend = nothing, color = ["grey" "grey83"], ylabel = [L"Certificate price [€/tCO$_2$]"], legend = (0.58,0.78), lw = 2, margin = 0mm, )  

        pprint = plot(ppp, size  = (500, 400),  xlim = [2023,(T_OPT_plot)],  titlefont = font(11),  title = Plot_title[C_nr], plot_titlefontsize = 12, topmargin = 2mm, gridlinewidth = 1, xtickfontsize=12,ytickfontsize=12, xguidefontsize=14,yguidefontsize=14,legendfontsize=12)
        

        filename_string  =  "C:/Users/VERBISTF/OneDrive - KU Leuven/PhD Flore/P4_negative_emissions/Figures/Plots/$(Case_name).svg"
        display(pprint)
        savefig(filename_string) 
        savefig(pprint, filename_string)


        return 
    end
end


function Joint_optimum_C1!()

    # NOTE 1: function works --> just make sure that you don't choose the discount rate too high, if the discountrate is too high , you have no banking taking place 
    # NOTE 2: also with removals there can be banking
    # NOTE 3: the EAC price is obtained by dividing the dual by the discount factor.

    model = Model(optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" =>  0))
    process_sets_parameters!(model, Parameters,Timeseries_data, T_OPT)


    model.ext[:variables] = Dict()
    model.ext[:constraints] = Dict()

    # Extract sets
    T = model.ext[:sets][:T] # set of years
    
    # Extract parameters
    DiscountFactor = model.ext[:parameters][:DiscountFactor] # Discount factor
    #e = model.ext[:parameters][:e]                           # Emissions emitted and not abated each year
    EAC_price = model.ext[:parameters][:EAC_price]           # EUA price 
    CDR_price = model.ext[:parameters][:CDR_price]           # CDR price



    ρ = model.ext[:parameters][:ρ] # Parameter controlling the price update step size in the ADMM algorithm
    CAP = model.ext[:parameters][:CAP] =  Dict(t =>  Timeseries_data[convert.(Int64,t+1), "ETS_cap [EUA/y]"] for t in T ) #Dict(t => BUDGET/(T_OPT+1) for t in T) # cap 
    N = model.ext[:parameters][:N] # number of participants in the ETS market 
    Budget = model.ext[:parameters][:Budget]

    # Define variables
    e = model.ext[:variables][:e] = @variable(model, e[T], lower_bound=0, upper_bound = EM_REF,  base_name="Emission allowances bought in year t by industrial fringe")
    r = model.ext[:variables][:r] = @variable(model, r[T], lower_bound=0,  upper_bound = EM_REF,  base_name="Emission allowances bought in year t by industrial fringe")

    ec = model.ext[:variables][:ec] = @variable(model, ec[T], lower_bound=0, base_name="Emission allowances bought in year t by industrial fringe")
    rc = model.ext[:variables][:rc] = @variable(model, rc[T], lower_bound=0, base_name="Emission allowances of which removal credits")

    # Constraints 
    Banking_e =  model.ext[:constraints][:Banking_e] = @constraint(model, [t=T],
    (sum(ec[i] for  i in 0:t)  >= sum(e[i]   for i in 0:t)))
    Banking_r =  model.ext[:constraints][:Banking_r] = @constraint(model, [t=T],
    (sum(r[i] for  i in 0:t)  >= sum(rc[i]   for i in 0:t)))
    Cap_con =  model.ext[:constraints][:Cap_con] = @constraint(model, [t=T],
    ec[t] == CAP[t] + rc[t] )

  
    No_removals =  model.ext[:constraints][:No_removals] = @constraint(model, [t=T],
    (rc[t] == 0))

    # Objective 
    
    model.ext[:objective] = @objective(model, Min, 
    # MACC
    sum((DiscountFactor[t])*(BETA*(EM_REF - e[t])^(2)/2  + (A_r*(r[t]/1) + B_r*(r[t])^(2)/EM_REF/2)) for t in T)) # NOTE: the /EM_REF is checked and correct
    
    
    optimize!(model);

    parameters_extract(model)
    variables_extract(model)

    
    EAC_surplus = [sum( ec[i]  - (e[i])  for i in 0:t) for t in model.ext[:sets][:T]]
    CDR_surplus = [sum((r[i] -rc[i])  for i in 0:t) for t in model.ext[:sets][:T]]
    EAC_price = [-dual.(model.ext[:constraints][:Cap_con][T_OPT]) for t in model.ext[:sets][:T]]

    Total_net_emissions = sum(e[t] - r[t] for t in model.ext[:sets][:T])
    Total_costs =sum((1/(1+DRc)^t)*(BETA*(EM_REF - e[t])^(2)/2  + (A_r*(r[t]/1) + B_r*(r[t])^(2)/EM_REF/2)) for t in  model.ext[:sets][:T])./Total_net_emissions
    # EA market plot^t
    Color = ["steelblue3" "teal"  "slategray3"]
    Alpha = [1 0.5 0.5]
    #EAC_price_2 = 1.5*Beta*(EM_REF .- EM_fringe).^Gamma
    p1 = stackedarea(Timeseries_data[1:1:(T_OPT+1),"Year"], [ec -rc], ylims = (0, 3000), labels = ["EAC bought by emitter"  "RC sold by remover"], title = "(a)", ylabel = "Number of certificates [Millions]", xlabel = "Years", legend = (0.32, -0.12), linecolor = Color, color = Color, alpha = Alpha, right_margin = 5mm)    
    plot!(twiny(p1), Timeseries_data[1:1:(T_OPT+1),"Year"], [Timeseries_data[1:1:(T_OPT+1), "ETS_cap [EUA/y]"]./1 EAC_surplus CDR_surplus], ylims = (0, 3000), label = ["EAC cap" "EAC surplus emitter" "RC surplus remover"], xformatter=_->"", yformatter=_->"", xlabel = "", ylabel="", legend = (0.02,0.98), color = ["black"  "navy" "black"], ls = [:solid :dash :dash], margin = 0mm)
    annotate!(2040, 100, Plots.text(L"%$(round(Total_costs_per_net, sigdigits = 3)) \textrm{€}/\textrm{tCO}_2-avg",  12))
    plot!(twinx(p1), Timeseries_data[1:1:(T_OPT+1),"Year"], [[EAC_price[t+1] for t in model.ext[:sets][:T]] [0/(1)^t  for t in model.ext[:sets][:T]]], ylims = (0, 450), yticks = 50:50:450, label = ["MAC (EAC price)" "MRC"], color = ["grey" "grey83"], ylabel = [L"Certificate price [€/tCO$_2$]"], legend = (0.58,0.98), lw = 2, margin = 0mm)


    Color = ["steelblue3" "slategray3"]
    Alpha = [1 0.5]

    p3 = stackedarea(Timeseries_data[1:1:(T_OPT+1),"Year"], [[e[i] for i in 0:1:T_OPT] [-r[i] for i in 0:1:T_OPT]], title = "(b)", labels = ["Emissions emitter" "Removals remover"], ylabel = "Emissions [Millions]", yguidefontsize= 10, xlabel = "Years", legend = (0.35, -0.12), linecolor = Color, color = Color, alpha = Alpha, left_margin = 5mm, grid= true, gridlinewidth = 0, margin = 0mm)
    # plot!(twiny(p3), xformatter=_->"", yformatter=_->"", xlabel = "", ylabel="", label = "", grid= false, gridlinewidth = 0, margin = 0mm)
    # plot!(twinx(p3), xformatter=_->"", yformatter=_->"", xlabel = "", ylabel="", label = "",grid= false, gridlinewidth = 0, margin = 0mm)
    annotate!(2040, 100, Plots.text(L"%$(round(Total_net_emissions./10^(3), sigdigits = 3))\textrm{GtCO}_2-tot",  12))


    l = @layout [grid(1,2)]
    p4 = plot(p1,p3, size  = (1000, 400), layout =l, titleloc = :right, titlefont = font(11), xlim = [2024,(2024+T_OPT)], plot_title = L"\parbox{50em}{\centering Joint Optimum \\}", plot_titlefontsize = 12, topmargin = 2mm)

    return model
end

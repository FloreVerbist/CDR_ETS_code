

@userplot StackedArea

# a simple "recipe" for Plots.jl to get stacked area plots
# usage: stackedarea(xvector, datamatrix, plotsoptions)
@recipe function f(pc::StackedArea)
    x, y = pc.args
    n = length(x)
    y = cumsum(y, dims=2)
    seriestype := :shape

    # create a filled polygon for each item
    for c=1:size(y,2)
        sx = vcat(x, reverse(x))
        sy = vcat(y[:,c], c==1 ? zeros(n) : reverse(y[:,c-1]))
        @series (sx, sy)
    end
end
using Plots; LaTeXStrings; 

pgfplotsx() # gr()

function variables_extract(model::Model)
    for key in keys(model.ext[:variables])
        V = value.(model.ext[:variables][key])
        @eval $key = $V
    end
    return
end

function parameters_extract(model::Model)
    for key in keys(model.ext[:parameters])
        P = model.ext[:parameters][key]
        @eval $key = $P
    end
end





function plot_pareto!(Dict_table, SOCIAL::Bool)
    ADMM = Dict()
    results = Dict()

    Dict_results_new = Dict()
    ETS_data = Dict() 
    compile_dictionary(Dict_results_new, ETS_data)
    if SOCIAL == true


        # for CO2_budget in Budget_range 
        #     REMOVALS = true
        #     model_joint = Joint_optimum_C0!(REMOVALS, CO2_budget)

        #     parameters_extract(model_joint)
        #     variables_extract(model_joint)

        #     Total_costs = sum(DiscountFactor[t]*(A_a*(1-e[t]./EM_REF) + B_a*((1-e[t]./EM_REF)).^(2)/2  + (A_r*(r[t]/EM_REF) + B_r*(r[t]/EM_REF).^(2)/2)) for t in 0:T_OPT).*EM_REF # NOTE: the /EM_REF is checked and correct

        #     Total_emissions  = sum(e[t] for t in 0:T_OPT) 

        #     push!(Dict_results_new[(CASE_NAME[1],"Total_net_emissions")], CO2_budget/1000)
        #     push!(Dict_results_new[(CASE_NAME[1],"Total_costs_per_net")], Total_costs./CO2_budget)
        #     push!(Dict_results_new[(CASE_NAME[1],"Total_emissions")], Total_emissions./1000)
            
        #     EM = [e[t] for t in 0:T_OPT]
        #     REM = [r[t] for t in 0:T_OPT]
        #     Color2 = ["snow3" "snow2"  "slategray3"]
        #     Alpha = [1 0.5 0.5]
        #     ppp = groupedbar(Timeseries_data[1:1:(T_OPT+1),"Year"], [EM -REM], bar_position =:stack, xticks = (2025:5:(2024+T_OPT)), ylims = (-600, 2500), yticks =-500:250:2500, labels = ["Emissions"  "Removals"], foreground_color_legend = nothing, title = "", ylabel = L"Emissions/Removals [MtCO$_2$/yr]", xlabel = "Years", legend = (0.02, 0.88), lw = 0.75, color = Color2, linecolor = [:gray31 :gray31], alpha = Alpha, right_margin = 5mm)
        #     plot!(twiny(ppp), Timeseries_data[1:1:(T_OPT+1),"Year"], [Timeseries_data[1:1:(T_OPT+1), "ETS_cap [EUA/y]"]./1], lw=[1 1], xticks = (2025:5:(2024+T_OPT)), ylims = (-600, 2500), label = ["ETS cap"],foreground_color_legend = nothing, xformatter=_->"", yformatter=_->"", xlabel = "", ylabel="", legend = (0.02,0.98), color = ["black"  "black" "black" "brown"], ls = [:dash :dash :dash :dash], margin = 0mm)
        #     plot!(twiny(ppp), Timeseries_data[1:1:(T_OPT+1),"Year"], [EM-REM], lw=0, xticks = (2025:5:(2024+T_OPT)), ylims = (-600, 2500), label = ["Net emissions"],foreground_color_legend = nothing, xformatter=_->"", yformatter=_->"", xlabel = "", ylabel="", legend = (0.02,0.93), color = ["black"  "black" "black" "brown"], ls = [:dash :dash :dash :dash], margin = 0mm, markershape = :xcross, markersize = 2)
    
        # end

    else

        for CO2_budget in Budget_range
            global Case_name = "C0-social_optimal"
            global C_nr = 1
            include("$(Case_name).jl")
            global RHO = RHO1 = RHO2 = RHO_vect[C_nr]
            global ADMM_max_it = 700

            EPSILON = 1
            define_results!(Parameters,results,ADMM)

            model_e = Model(optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" =>  0))
            model_r = Model(optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" =>  0))
            model_ccb = Model(optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" =>  0))


            # Step 3b: process parameters
            global BUDGET = CO2_budget
            print(BUDGET)
            process_sets_parameters!(model_e, Parameters,Timeseries_data, T_OPT)
            process_sets_parameters!(model_r, Parameters,Timeseries_data, T_OPT)
            process_sets_parameters!(model_ccb, Parameters,Timeseries_data, T_OPT)


            # calculate equilibrium 
            ADMM["walltime"] = @elapsed Dict_results_new = ADMM_price_based_m!(Case_name,results, ADMM, model_e, model_r, model_ccb, Dict_results_new, LINEAR)  
        
        
        end
        CSV.write("Pareto_line.csv", Dict_results_new)
    end
        Dict_results_new = CSV.read("Pareto_line.csv", DataFrame)

        xmin = 8
        xmax = 40
        ymin = 50
        ymax = 750
        xstep = 4
        ystep = 50
        
        # xmin = 5
        # xmax = 25
        # ymin = 650
        # ymax = 1650
        # xstep = 5
        # ystep = 100
        p1= plot(xlims= (xmin, xmax), ylims=(ymin, ymax), foreground_color_legend = nothing, legend = (0.35,0.97), xtickfontsize=12,ytickfontsize=12, xguidefontsize=14,yguidefontsize=14)
        MS = [:star6, :rect, :circle, :pentagon,  :dtriangle, :utriangle, :diamond, :rtriangle, :ltriangle]
        MC = ["springgreen4", "black", "steelblue4", "deepskyblue3", "slateblue3", "thistle", "mediumorchid2", "chocolate1", "brown3"]
        MSW = [3, 1, 1, 1, 1, 1, 1, 1,1]
        MSC = ["springgreen4", "black", "black", "black", "black", "black", "black", "black", "black"]

        #MS = [:circle, :circle, :circle, :circle, :circle, :rtriangle]
        for (i, Plot_title_name) in enumerate(Plot_title_ordered[1:end])
            if i != 100
            C_nr = i
            p1 = scatter!([Dict_table[Dict_table[:, "Case_name"] .== Plot_title_name,"Total_net_emissions"]], [Dict_table[Dict_table[:, "Case_name"] .== Plot_title_name, "Total_costs_per_net"][end]], labels= Plot_title_ordered[C_nr], ylabel = L"Total cost per net emission [€/tCO$_2$]", xlabel=L"Carbon budget (total net emissions) [GtCO$_2$]", legendfonthalign = :left, markershape = MS[C_nr], markercolor = MC[C_nr],  markerstrokewidth= MSW[C_nr], markerstrokecolor = MSC[C_nr], markersize = 6)#,  legend = (0.9, 1))
            else 
            end
        end
        social_pareto_front_string = Dict_results_new[Dict_results_new[!,1] .== "(\"C0-social_optimal\", \"Total_costs_per_net\")", 2][1][1:end]
        social_pareto_front_vector = eval(Meta.parse(replace(social_pareto_front_string, "Any" => "")))
        carbon_budget_front_string = Dict_results_new[Dict_results_new[!,1] .== "(\"C0-social_optimal\", \"Total_net_emissions\")", 2][1][1:end]
        carbon_budget_front_vector = eval(Meta.parse(replace(carbon_budget_front_string, "Any" => "")))
        p1 = plot!([carbon_budget_front_vector], [social_pareto_front_vector], labels = "Pareto Front", legendfonthalign = :left, color = "black", ls="dash", lw=1)
        # plot!([Dict_results_new[(CASE_NAME[1], "Total_net_emissions")][1:end]], [Dict_results_new_initial[(CASE_NAME[1], "Total_costs_per_net")][1:end]], labels = "Pareto Front", legendfonthalign = :left, color = "blue", ls="dash", lw=1)
        plot!(twinx(p1), [0 0], [0 0], xlims= (xmin, xmax), ylims=(ymin, ymax), xticks = xmin:xstep:xmax, yticks = ymin:ystep:ymax,label="",  color = :black, xformatter=_->"", yformatter=_->"") 
        plot!(twiny(p1), [0 0], [0 0], xlims= (xmin, xmax), ylims=(ymin, ymax),  xticks = xmin:xstep:xmax, yticks = ymin:ystep:ymax, label="", color = :black, xformatter=_->"", yformatter=_->"")
        plot!(grid = true,  size  = (600, 500), gridalpha=0.03, gridlinewidth=1,xlims= (xmin, xmax), ylims=(ymin, ymax),  xticks = xmin:xstep:xmax, yticks = ymin:ystep:ymax)

        filename_string  =  "C:/Users/VERBISTF/OneDrive - KU Leuven/PhD Flore/P4_negative_emissions/Figures/Plots/pareto_plot.svg"
        savefig(filename_string) 
        savefig(p1, filename_string)

    return display(p1)
end

function ccb_budget_plot(Dict_table)
    # average 
    MS = [:star6, :rect, :circle, :pentagon,  :dtriangle, :utriangle, :diamond, :rtriangle, :ltriangle]
    MC = ["springgreen4", "black", "steelblue4", "deepskyblue3", "slateblue3", "thistle", "mediumorchid2", "chocolate1", "brown3"]
    MSW = [3, 1, 1, 1, 1, 1, 1, 1,1]
    MSC = ["springgreen4", "black", "black", "black", "black", "black", "black", "black", "black"]
    X_labels = [L"\parbox{7em}{\centering No \\ Removals}", L"\parbox{7em}{\centering Target-based}", L"\parbox{7em}{\centering Net Emission \\ Cap}",    L"\parbox{7em}{\centering Removal \\ Obligation}", L"\parbox{7em}{\centering Equivalence \\ Obligation}",  L"\parbox{7em}{\centering Conditional \\ Integration}",  L"\parbox{7em}{\centering Gross Emission \\ Cap}",   L"\parbox{7em}{\centering Gross cap \\ Obligation}"]
    CCB_budget = Dict_table[:, "Total_budget"]
    textcolor = [:black :black :black :black :black :black :black :black :black]
    
    legend_p = (0.02, 0.98)

    pp = bar([1:8], [CCB_budget[2:end]], label = nothing, legend = nothing, color = nothing, linecolor = nothing, ylim = (-50, 260), xlim = (8.5,20.75), xrotation=90, ylabel = L"Profits or losses CCB per net emission [€/tCO$_2$]") #, color = MC)
    plot!(twinx(pp),  xformatter=_->"", yformatter=_->"", xlabel = "", ylabel="", ylim = (-50, 260), xlim = (8.5,20.75))
    plot!(twiny(pp),  xformatter=_->"", yformatter=_->"", xlabel = "", ylabel="", ylim = (-50, 260), xlim = (8.5,20.75), grid= :y, gridlinewidth = 1.0 )
    for i in 1:length(X_labels)
        bar!([X_labels[i]], [CCB_budget[2:end][i]], color=MC[2:end][i], linecolor = MSC[2:end][i], alpha = [0.9],  label = nothing, xrotation=90, legend = nothing)
        annotate!(X_labels[i], CCB_budget[2:end][i], text(round(CCB_budget[2:end][i], sigdigits =4), textcolor[i], 8, :bottom))
    end
    
    pprint = plot(pp, size  = (500, 400), title = "CCB balance", yticks = -50:50:250, titlefont = font(11), xlim = (8.5,20.75),  plot_titlefontsize = 6, topmargin = 2mm, xtickfontsize = 10)

    filename_string  =  "C:/Users/VERBISTF/OneDrive - KU Leuven/PhD Flore/P4_negative_emissions/Figures/Plots/revenues_CCB.svg"
    display(pprint)
    savefig(filename_string) 
    savefig(pprint, filename_string)


    return display(pprint)
end

function carbon_lock_in_plot(Dict_table)
    # average spending per gross emissions
    #Dict_table = CSV.read("Table_for_paper_sept.csv", DataFrame)
    ymax = 60
    ymin = - 40
    MS = [:star6, :rect, :circle, :pentagon,  :dtriangle, :utriangle, :diamond, :rtriangle, :ltriangle]
    MC = ["springgreen4", "black", "steelblue4", "deepskyblue3", "slateblue3", "thistle", "mediumorchid2", "chocolate1", "brown3"]
    MSW = [3, 1, 1, 1, 1, 1, 1, 1,1]
    MSC = ["springgreen4", "black", "black", "black", "black", "black", "black", "black", "black"]
    MSC_net = ["springgreen4", "black", "black", "black", "black", "black", "black", "black", "black"]

    X_labels = [L"\parbox{7em}{\centering No \\ Removals}", L"\parbox{7em}{\centering Target-based}", L"\parbox{7em}{\centering Net Emission \\ Cap}",    L"\parbox{7em}{\centering Removal \\ Obligation}", L"\parbox{7em}{\centering Equivalence \\ Obligation}",  L"\parbox{7em}{\centering Conditional \\ Integration}",  L"\parbox{7em}{\centering Gross Emission \\ Cap}",   L"\parbox{7em}{\centering Gross cap \\ Obligation}"]
    Emissions = Dict_table[:, "Total_emissions"]
    Removals = Dict_table[:, "Total_emissions"] - Dict_table[:, "Total_net_emissions"]
    Net =  Dict_table[:, "Total_net_emissions"]
    textcolor = [:black :black :black :black :black :black :black :black :black]
    
    legend_p = (0.02, 0.98)

    pp = bar([1:8], [Emissions[2:end] -Removals[2:end]], label = nothing, legend = nothing, color = nothing, linecolor = nothing, ylim = (ymin, ymax), xlim = (8.5,20.75), xrotation=90, ylabel = L"Total emissions/removals [GtCO$_2$]") #, color = MC)
    plot!(twinx(pp),  xformatter=_->"", yformatter=_->"", xlabel = "", ylabel="", ylim = (ymin, ymax), xlim = (8.5,20.75))
    plot!(twiny(pp),  xformatter=_->"", yformatter=_->"", xlabel = "", ylabel="", ylim = (ymin, ymax), xlim = (8.5,20.75), grid= :y, gridlinewidth = 1.0 )
    for i in 1:length(X_labels)
        bar!([X_labels[i]], [Emissions[2:end][i] -Removals[2:end][i]], alpha = [0.9 0.4], color=MC[2:end][i], linecolor = MSC[2:end][i],  label = nothing, xrotation=90, legend = nothing)
        scatter!([X_labels[i]], [Net[2:end][i]], markershape = :hline, markercolor = MSC_net, linestyle = :dash )
        annotate!(X_labels[i], Emissions[2:end][i], text(round(Emissions[2:end][i], sigdigits =3), textcolor[i], 10, :bottom))
        annotate!(X_labels[i], Net[2:end][i], text(round(Net[2:end][i], sigdigits =3), textcolor[i], 10, :bottom))
        annotate!(X_labels[i], -Removals[2:end][i], text(round(Removals[2:end][i], sigdigits =3), textcolor[i], 10, :bottom))

    end
    annotate!(10.7, 58, text("Total emissions", :black, 10, :top))
    annotate!(10.7, 53, text("Total removals", :black, 10, :top))
    annotate!(10.65, 48, text("Net emissions", :black, 10, :top))

    annotate!(9, 47, text("_", :black, 8, :top))
    scatter!([9], [55], markershape = :rect, alpha = [0.9], color = [:grey])
    scatter!([9], [50], markershape = :rect, alpha = [0.4], color = [:grey])
    # plot!( [(9,375), (9,400)], arrow = arrow(4, 4), color =:black)
   
    pprint = plot(pp, size  = (500, 400), title = "Emission balance", yticks = ymin:10:ymax, titlefont = font(11), xlim = (8.5,20.75),  plot_titlefontsize = 6, topmargin = 2mm, xtickfontsize = 10)

    filename_string  =  "C:/Users/VERBISTF/OneDrive - KU Leuven/PhD Flore/P4_negative_emissions/Figures/Plots/emission_balance.svg"
    display(pprint)
    savefig(filename_string) 
    savefig(pprint, filename_string)


    return display(pprint)
end

function emitter_budget_plot(Dict_table)
    # average spending per gross emissions

    MS = [:star6, :rect, :circle, :pentagon,  :dtriangle, :utriangle, :diamond, :rtriangle, :ltriangle]
    MC = ["springgreen4", "black", "steelblue4", "deepskyblue3", "slateblue3", "thistle", "mediumorchid2", "chocolate1", "brown3"]
    MSW = [3, 1, 1, 1, 1, 1, 1, 1,1]
    MSC = ["springgreen4", "black", "black", "black", "black", "black", "black", "black", "black"]
    X_labels = [L"\parbox{7em}{\centering No \\ Removals}", L"\parbox{7em}{\centering Target-based}", L"\parbox{7em}{\centering Net Emission \\ Cap}",    L"\parbox{7em}{\centering Removal \\ Obligation}", L"\parbox{7em}{\centering Equivalence \\ Obligation}",  L"\parbox{7em}{\centering Conditional \\ Integration}",  L"\parbox{7em}{\centering Gross Emission \\ Cap}",   L"\parbox{7em}{\centering Gross cap \\ Obligation}"]
    Emitter_budget = Dict_table[:, "Emitter_cost_avg"]
    textcolor = [:black :black :black :black :black :black :black :black :black]
    
    legend_p = (0.02, 0.98)

    pp = bar([1:8], [Emitter_budget[2:end]], label = nothing, legend = nothing, color = nothing, linecolor = nothing, ylim = (-50, 450), xlim = (8.5,20.75), xrotation=90, ylabel = L"Costs emitters per gross emissions [€/tCO$_2$]") #, color = MC)
    plot!(twinx(pp),  xformatter=_->"", yformatter=_->"", xlabel = "", ylabel="", ylim = (-50, 450), xlim = (8.5,20.75))
    plot!(twiny(pp),  xformatter=_->"", yformatter=_->"", xlabel = "", ylabel="", ylim = (-50, 450), xlim = (8.5,20.75), grid= :y, gridlinewidth = 1.0 )
    for i in 1:length(X_labels)
        bar!([X_labels[i]], [Emitter_budget[2:end][i]], color=MC[2:end][i], linecolor = MSC[2:end][i], alpha = [0.9],  label = nothing, xrotation=90, legend = nothing)
        annotate!(X_labels[i], Emitter_budget[2:end][i], text(round(Emitter_budget[2:end][i], sigdigits =4), textcolor[i], 8, :bottom))
    end
    annotate!(8.75, 350, text(round(Emitter_budget[2:end][1], sigdigits =4), :black, 8, rotation = 90, :top))
    plot!( [(9,375), (9,400)], arrow = arrow(4, 4), color =:black)
    #plot!([8.3, 17.5], [Dict_table[2:end, "Total_budget"][1], Dict_table[2:end, "Total_budget"][1]], color = :black, linestyle = :dash,  ylim = (-50, 450), legend = nothing, label = nothing)
    # bar!([X_labels[1]], [0], color=MC[2], alpha = 0.9, label = "Reference" , legend = legend_p, legendfonthalign = :left)
    # bar!([X_labels[2]], [0], color=MC[3], alpha = 0.9, label = "Low ETS price & high RC price" , legend = legend_p)
    # bar!([X_labels[3]], [0], color=MC[4], alpha = 0.9, label = "Low ETS price" , legend = legend_p)
    # bar!([X_labels[4]], [0], color=MC[5], alpha = 0.9, label = "High ETS price & high number EAC" , legend = legend_p)
    # bar!([X_labels[6]], [0], color=MC[7], alpha = 0.9, label = "High ETS price & low net emissions" , legend = legend_p)
    pprint = plot(pp, size  = (500, 400), title = "Emitter balance", yticks = -50:50:450, titlefont = font(11), xlim = (8.5,20.75),  plot_titlefontsize = 6, topmargin = 2mm, xtickfontsize = 10)

    filename_string  =  "C:/Users/VERBISTF/OneDrive - KU Leuven/PhD Flore/P4_negative_emissions/Figures/Plots/spendings_emitters.svg"
    display(pprint)
    savefig(filename_string) 
    savefig(pprint, filename_string)


    return display(pprint)
end

# function ccb_emitter_balance(Dict_table)

#     MS = [:star6, :rect, :circle, :pentagon,  :dtriangle, :utriangle, :diamond, :rtriangle, :ltriangle]
#     MC = ["springgreen4", "black", "steelblue4", "deepskyblue3", "slateblue3", "thistle", "mediumorchid2", "chocolate1", "brown3"]
#     MSW = [3, 1, 1, 1, 1, 1, 1, 1,1]
#     MSC = ["springgreen4", "black", "black", "black", "black", "black", "black", "black", "black"]
#     X_labels = [L"\parbox{7em}{\centering No \\ Removals}", L"\parbox{7em}{\centering Target-based}", L"\parbox{7em}{\centering Net Emission \\ Cap}",    L"\parbox{7em}{\centering Removal \\ Obligation}", L"\parbox{7em}{\centering Equivalence \\ Obligation}",  L"\parbox{7em}{\centering Conditional \\ Integration}",  L"\parbox{7em}{\centering Gross Emission \\ Cap}",   L"\parbox{7em}{\centering Gross cap \\ Obligation}"]
#     Emitter_budget = Dict_table[2:end, "Emitter_cost_avg"]
#     textcolor = [:black :black :black :black :black :black :black :black :black]
#     CCB_budget = Dict_table[2:end, "Total_budget"]


    
#     # Creating the x-axis points (as you indicated, from 9 to 16 to match your xlim)
#     x = 9:16
    
#     # Create grouped bar plot using the correct syntax and settings
#     pp = groupedbar(x, [Emitter_budget CCB_budget], 
#         label = ["Emitter Budget" "CCB Budget"],   # You can remove this if you want default colors
#         ylim = (-50, 450), 
#         xlim = (8.5, 16.5), 
#         xrotation = 90, 
#         ylabel = L"Costs emitters per gross emissions [€/tCO_2]", 
#         xticks = (x, X_labels), 
#         color = MC[2:end]
#     )
    
#     return display(pprint)

# end

function learning_plot(Dict_table)
    Dict_table = CSV.read("Table_for_paper_sept.csv", DataFrame)
    YMAX = 250
    XMAX = 50
    pp = plot(title="", legend = (0.03, 0.98), foreground_color_legend = nothing, gridlinewidth = 2, xtickfontsize=12,ytickfontsize=12, xguidefontsize=14,yguidefontsize=14)

    MS = [:star6, :rect, :pentagon, :diamond,  :circle, :utriangle, :dtriangle, :rtriangle, :ltriangle]
    MC = ["springgreen4", "black", "steelblue4", "deepskyblue3", "slateblue3", "thistle", "mediumorchid2", "chocolate1", "brown3"]
    MC_sense = ["springgreen4", "black", "steelblue4", "deepskyblue3", "slateblue3", "thistle", "mediumorchid2", "chocolate1", "brown3"]
    MSW = [3, 2, 2, 2, 2, 2, 2, 2, 2]
    MSC = ["springgreen4", "black", "black", "black", "black", "black", "black", "black", "black"]
    ordered_Cnr_pareto = [1, 2, 3, 4, 7, 8, 6, 5, 9]

    SENSITIVITY = true
    if SENSITIVITY == true

        Sensitivity_tab = CSV.read("Sensitivity_dict_table.csv", DataFrame)
        Sensitivity_tab[!, 2] = [eval(Meta.parse(replace(Sensitivity_tab[i,2], "Any" => ""))) for i in 1:length(Sensitivity_tab[!,1])] 
        data = Sensitivity_tab 

        for (i, Case_name) in enumerate(CASE_NAME_ordered[3:end])
            C_nr = i+2
            Total_gross = data[data[:,1] .== "(\"$(CASE_NAME_ordered[C_nr])\", \"Total_emissions\")",2] 
            Total_net = data[data[:,1] .== "(\"$(CASE_NAME_ordered[C_nr])\", \"Total_net_emissions\")",2] 
            Fraction_reality = [EQ_range[m] for m in Equivalence_vector]
            L_R_gains = data[data[:,1] .== "(\"$(CASE_NAME_ordered[C_nr])\", \"Learning_gains_per_rem\")",2] 
            Total_gross_removals = ( vcat(Total_gross...) .- vcat(Total_net...))./Fraction_reality
            scatter!([Total_gross_removals], [L_R_gains],  markershape = MS[C_nr], markercolor = MC_sense[C_nr], markeralpha = 0.3, markersize = 4, markerstrokewidth = 0, markerstrokecolor = MSC[C_nr], xlims = (-1, XMAX), xticks =0:5:XMAX,ylim = (-5, YMAX), yticks = 0:25:YMAX, label="")#,  legend = (0.9, 1))
        
        end
    else
    end

    scatter!([17.5], [65], color = :white, markeralpha = 0.6,  markerstrokecolor = :white, markersize = 60, label = "")

    # PLotting markers
    for (i, Case_name) in enumerate(CASE_NAME[1:end])
        C_nr = i
        scatter!([( Dict_table[C_nr,"Total_emissions"] - Dict_table[C_nr,"Total_net_emissions"])./FRAC_R_vec[ordered_Cnr_pareto[C_nr]]], [Dict_table[C_nr, "Learning_gains_per_rem"][end]], labels= Plot_title_ordered[C_nr], ylabel = L"Learning gains per gross removals [€/tCO$_2$]", xlabel=L"Cumulative gross removals [GtCO$_2$]", legendfonthalign = :left, background_color_legend = nothing, markershape = MS[C_nr], markercolor = MC[C_nr], markersize = 7, markerstrokewidth = MSW[C_nr], markerstrokecolor = MSC[C_nr], xlims = (-1, XMAX), xticks =0:5:XMAX,ylim = (-5, YMAX), yticks = 0:25:YMAX)#,  legend = (0.9, 1))
    end
    # Plotting lines: 
    for (i, Case_name) in enumerate(CASE_NAME[1:end])
        C_nr = i
        plot!(twinx(pp), [-5, ( Dict_table[C_nr,"Total_emissions"] - Dict_table[C_nr,"Total_net_emissions"])./FRAC_R_vec[ordered_Cnr_pareto[C_nr]]], [Dict_table[C_nr, "Learning_gains_per_rem"][end], Dict_table[C_nr, "Learning_gains_per_rem"][end]], label="", ylabel="",  xformatter=_->"", yformatter=_->"", xlabel = "", alpha = 0.35, color = "black", lw = 1, ls = :dash,  xlims = (-1, XMAX), xticks =0:5:XMAX, ylim = (-5, YMAX))
        plot!(twinx(pp), [( Dict_table[C_nr,"Total_emissions"] - Dict_table[C_nr,"Total_net_emissions"])./FRAC_R_vec[ordered_Cnr_pareto[C_nr]], ( Dict_table[C_nr,"Total_emissions"] - Dict_table[C_nr,"Total_net_emissions"])./FRAC_R_vec[ordered_Cnr_pareto[C_nr]]], [-5, Dict_table[C_nr, "Learning_gains_per_rem"][end]], alpha = 0.35, label="", ylabel="",  xformatter=_->"", yformatter=_->"", xlabel = "", color = "black", lw = 1, ls = :dash, xlims = (-1, XMAX), xticks =0:5:XMAX, ylim = (-5, YMAX))
    end


    #Not used:
    plot!(twiny(pp), sort(( Dict_table[2:end,"Total_emissions"] - Dict_table[2:end,"Total_net_emissions"])./FRAC_R_vec[ordered_Cnr_pareto[C_nr]]; alg=QuickSort),sort(Dict_table[2:end, "Learning_gains_per_rem"]; alg=QuickSort), label="", ylabel="",  xformatter=_->"", yformatter=_->"", xlabel = "", color = "mediumseagreen", lw = 2, alpha = 0.0, xlims = (-1, XMAX), xticks =0:5:XMAX, ylim = (-5, YMAX))

    plot!(twinx(pp),  xformatter=_->"", yformatter=_->"", xlabel = "", ylabel="", xlims = (-1, XMAX), ylim = (-5, YMAX), xticks =0:5:XMAX, yticks = 0:25:YMAX, grid = true, gridalpha=0.05, gridlinewidth=1)


    filename_string  =  "C:/Users/VERBISTF/OneDrive - KU Leuven/PhD Flore/P4_negative_emissions/Figures/Plots/learning_removal_plot.svg"
    savefig(filename_string) 
    savefig(pp, filename_string)


    return display(pp)
end


function case_data!(Dict1::Dict, Dict2::Dict, Dict3::Dict, Dict4::Dict)
    variable_names = ["Total_economic_surplus", "Total_costs", "Total_budget", "Total_net_emissions", "Total_emissions", "ETS_max", "ETS_var"]
    data = [[Dict1[var] for var in variable_names], [Dict2[var] for var in variable_names], [Dict3[var] for var in variable_names], [Dict4[var] for var in variable_names]]
    return data
end



function sensitivity_analysis()

    EPSILON = 5.0
    M_r_range = 0.15:0.15:0.45
    D_r_range = [2000, 2500, 3000]
    D_a_range = [50, 100, 200]
    EQ_range = [1, 0.85]
    scenario_level = ["L" "M" "H"]
    scenario_level_2 = ["yes" "no"]
    scenario_vector = []
    D_r_vector = []
    D_a_vector = []
    M_r_vector = []
    Equivalence_vector = []
    ADMM_max_it = 800


    Dict_results_SA = Dict()
    ETS_data_SA = Dict() 
    compile_dictionary(Dict_results_SA, ETS_data_SA)

    for M_R_level in 1:length(M_r_range)
        for D_R_level in 1:length(D_r_range)
            for D_A_level in 1:length(D_a_range)
                for EQ_level in 1:length(EQ_range)
                   
                    Scenario = "Learning R= $(scenario_level[M_R_level]), Cost R = $(scenario_level[D_R_level]), Cost A = $(scenario_level[D_A_level]), Equivalence = $(scenario_level_2[EQ_level])"
                    push!(scenario_vector, Scenario) #, "Learning A= $(scenario_level[M_A])")
                    push!(D_r_vector, D_R_level) # MISTAKE HERE!!!! --> corrected in paper
                    push!(D_a_vector, D_A_level)
                    push!(M_r_vector, M_R_level)
                    push!(Equivalence_vector, EQ_level)
                    

                    # for (index, Case_name) in enumerate(CASE_NAME[3:end])
                    #     C_nr = index + 2
                        ADMM_max_it = 800
                        # try 

                        Case_name = CASE_NAME[C_nr]
                        RHO = RHO1 = RHO2 = RHO_vect[C_nr]
                        ###--------------------------------------------------------------
                        results_SA = Dict()
                        ADMM_SA = Dict()
                        global BUDGET = Carbon_budget
                    
                        print("$(Case_name)  --  ")
                        include("$(Case_name).jl")
                    
                        define_results!(Parameters,results_SA,ADMM_SA)
                    
                        global model_e = Model(optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" =>  0))
                        global model_r = Model(optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" =>  0))
                        global model_ccb = Model(optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" =>  0))
                    
                        global M_R = M_r_range[M_R_level]
                        global D_R = D_r_range[D_R_level]
                        global D_A = D_a_range[D_A_level]
                        global FRAC_R = EQ_range[EQ_level]
                        # Step 3b: process parameters
                        process_sets_parameters!(model_e, Parameters,Timeseries_data, T_OPT)
                        process_sets_parameters!(model_r, Parameters,Timeseries_data, T_OPT)
                        process_sets_parameters!(model_ccb, Parameters,Timeseries_data, T_OPT)
                    
                        # recalibration parameters 

                    
                        # calculate equilibrium 
                        ADMM_SA["walltime"] = @elapsed Dict_results_SA = ADMM_price_based_m!(Case_name,results_SA, ADMM_SA, model_e, model_r, model_ccb, Dict_results_SA, ETS_data_SA, LINEAR)  
                        #print("Solved for $(Scenario) and case $(Case_name)") 
                        # catch
                        #     print("Error for $(Scenario) and case $(Case_name)") 
                        # end
                    
                    # end 
                    

                end
            end
        end
    end

    Dict_results_SA_cond_2 = Dict_results_SA
    ETS_data_SA_cond_2 = ETS_data_SA

    CSV.write("Sensitivity_dict_table_cond_2.csv", Dict_results_SA_cond_2)
    CSV.write("Sensitivity_table_ETS_cond_2.csv", ETS_data_SA_cond_2)


    # CSV.write("Sensitivity_dict_table.csv", Dict_results_SA)
    # CSV.write("Sensitivity_table_ETS.csv", ETS_data_SA)
   
   
    # CSV.write("Scenario_vector.csv", scenario_vector)
end

function postprocessing()
    # To involve the no removals: 
    #########################################################""
    data_original = CSV.read("C:/Users/VERBISTF/OneDrive - KU Leuven/PhD Flore/P4_negative_emissions/Scripts/Script_P2_CDR_ETS/Script_P2_CDR_ETS/SAVED_outputs/Sensitivity_dict_table.csv", DataFrame)
    data_no_removals = CSV.read("Sensitivity_dict_table_no_removals.csv", DataFrame)
    for (i, value) in enumerate(Table_column_names)
        #data_original[data_original[:,1] .== "(\"$(CASE_NAME[2])\", \"$(Table_column_names[i])\")",2] =          data_no_removals[data_no_removals[:,1] .== "(\"$(CASE_NAME[2])\", \"$(Table_column_names[i])\")",2]
        data_original[data_original[:,1] .== "(\"$(CASE_NAME[2])\", \"$(Table_column_names[i])\")",2] =          data_no_removals[data_no_removals[:,1] .== "(\"$(CASE_NAME[2])\", \"$(Table_column_names[i])\")",2]

    end

    # CSV.write("C:/Users/VERBISTF/OneDrive - KU Leuven/PhD Flore/P4_negative_emissions/Scripts/Script_P2_CDR_ETS/Script_P2_CDR_ETS/SAVED_outputs/Sensitivity_dict_table.csv", data_original)
    #####################################################

    M_r_range = 0.15:0.15:0.45
    D_r_range = [2000, 2500, 3000]
    D_a_range = [50, 100, 200]
    EQ_range = [1, 0.85]
    scenario_level = ["L" "M" "H"]
    scenario_level_2 = ["yes" "no"]
    scenario_vector = []
    D_r_vector = []
    D_a_vector = []
    M_r_vector = []
    Equivalence_vector = []
    
    for M_R_level in 1:length(M_r_range)
        for D_R_level in 1:length(D_r_range)
            for D_A_level in 1:length(D_a_range)
                for EQ_level in 1:length(EQ_range)
                   
                    Scenario = "Learning R= $(scenario_level[M_R_level]), Cost R = $(scenario_level[D_R_level]), Cost A = $(scenario_level[D_A_level]), Equivalence = $(scenario_level_2[EQ_level])"
                    push!(scenario_vector, Scenario) #, "Learning A= $(scenario_level[M_A])")
                    push!(D_r_vector, D_R_level) # MISTAKE HERE --> corrected in paper
                    push!(D_a_vector, D_A_level)
                    push!(M_r_vector, M_R_level)
                    push!(Equivalence_vector, EQ_level)
                    
                end
            end
        end
    end

end

function scenario_bars!(data)
    Sensitivity_tab = CSV.read("Sensitivity_dict_table.csv", DataFrame)
    Sensitivity_tab = data_original # This contains also the no removals data - as it was added later --> see processing function
    Sensitivity_tab[!, 2] = [eval(Meta.parse(replace(Sensitivity_tab[i,2], "Any" => ""))) for i in 1:length(Sensitivity_tab[!,1])] 
    data = Sensitivity_tab 
    Case_data = CSV.read("Table_for_paper_sept.csv", DataFrame)
    #Case_data[:, "Total_budget_per_net"] = Case_data[:, "Total_budget"]
    variable_names = ["Total_net_emissions", "Total_emissions", "Emitter_cost_avg", "Total_budget_per_net",  "Learning_gains_per_rem", "Total_costs_per_net"]
    plot_title = ["Total net emissions", "Total gross emissions", "Emitter balance", "CCB balance", "Learning gains removals", "Total cost per net emission"]
    plot_title = ["","","","","","","","",""]
    y_title = [L"Net emissions [GtCO$_2$]",  L"Gross emissions [GtCO$_2$]", L"Emitters' costs per gross emission [€$$/tCO$_2$]", L"Income per net emission [€$$/tCO$_2$]", L"Economic gains per gross removal [€$$/tCO$_2$]",  L"Total costs per net emission [€$$/tCO$_2$]"]
    Model_names = [L"\parbox{7em}{\centering Social \\ Optimum}", L"\parbox{7em}{\centering No \\ Removals}", L"\parbox{7em}{\centering Target-based}", L"\parbox{7em}{\centering Net Emission \\ Cap}",    L"\parbox{7em}{\centering Removal \\ Obligation}", L"\parbox{7em}{\centering Equivalence \\ Obligation}",  L"\parbox{7em}{\centering Conditional \\ Integration}",  L"\parbox{7em}{\centering Gross Emission \\ Cap}",   L"\parbox{7em}{\centering Gross cap \\ Obligation}"]
    MC = ["springgreen4", "black", "steelblue4", "deepskyblue3", "slateblue3", "thistle", "mediumorchid2", "chocolate1", "brown3"]
    ordered_Cnr = [2, 3, 4, 7, 8, 6, 5, 9]
    scattercolors = ["black" "gray73" "white"] #range(colorant"red", stop=colorant"green", length=length(data[1][1]))
    scatteralpha = [0.9 0.5 0.1]
    scattermarkershape = [:circle, :rect, :utriangle]
    scattersize = [2,3,4]
    markerwidth = [0, 1]
    ymin_vec = [10, 10, 100, -500, -20, 0]
    ymax_vec = [35, 60, 800, 800, 250, 700]
    # index_nr = 6
    # i = variable_names[index_nr]
    for (index_nr,i) in enumerate(variable_names)

        ymax = maximum(hcat([data[data[!,1] .== "(\"$(CASE_NAME[n])\", \"$(i)\")", 2][1][1:end] for n in ordered_Cnr ]...))
        ymin = minimum(hcat([data[data[!,1] .== "(\"$(CASE_NAME[n])\", \"$(i)\")", 2][1][1:end] for n in ordered_Cnr ]...))
        ydelta = 0.05 * (ymax - ymin)
        YMAX = ymax + ydelta
        YMIN = ymin - ydelta

        p = boxplot(reshape(hcat([data[data[!,1] .== "(\"$(CASE_NAME[n])\", \"$(i)\")", 2][1][1:end] for n in ordered_Cnr ]...), 54,length(ordered_Cnr)), xticks = ([m for m in 1:length(ordered_Cnr)], [Model_names[n] for n in 2:(1+length(ordered_Cnr))]), alpha = 0.5, label = "", title = "$(plot_title[index_nr])", ylabel = y_title[index_nr], xlabel = "", legend = false, xrotation=90, color = hcat([MC[m] for m in 2:(1+length(ordered_Cnr))]...), ylim = (YMIN, YMAX), xtickfontsize=12,ytickfontsize=12, xguidefontsize=14,yguidefontsize=14)

        for j in 2:1:length(data[!,2][1])-1
            scatter!([m for m in 1:length(ordered_Cnr)],[data[data[!,1] .== "(\"$(CASE_NAME[n])\", \"$(i)\")", 2][1][j] for n in ordered_Cnr ],  color = [scattercolors[M_r_vector[j]] for m in 1:length(ordered_Cnr)], markershape = [scattermarkershape[D_r_vector[j]] for m in 1:length(ordered_Cnr)], markersize = [scattersize[D_a_vector[j]] for m in 1:length(ordered_Cnr)], markerstrokewidth =[markerwidth[Equivalence_vector[j]] for m in 1:length(ordered_Cnr)], ylim = (YMIN, YMAX))
        end
        scatter!([m for m in 1:length(ordered_Cnr)], [Case_data[(m+1), i] for m in 1:length(ordered_Cnr)], markershape = [:xcross], markeralpha = [0.5], markersize = [8], markercolor = [:black], markerstrokewidth =[3])
        plot!(twinx(p), label="",  color = :white, xformatter=_->"", yformatter=_->"", ylim = (YMIN, YMAX)) 
        plot!(twiny(p), label="",  color = :white, xformatter=_->"", yformatter=_->"", gridlinewidth = 0, grid =:hide, gridalpha = 0, gridstyle = :dash, ylim = (YMIN, YMAX)) 
        plot!([2.5, 2.5], [ymin_vec[index_nr], ymax_vec[index_nr]], ls = :dash, color = :black, label = "")
        plot!([6.5, 6.5], [ymin_vec[index_nr], ymax_vec[index_nr]], ls = :dash, color = :black, label = "")
        plot!(grid = true, size  = (600, 400), ylim = (ymin_vec[index_nr], ymax_vec[index_nr])) # ylim = (0, 700)
        display(p)
        # filename_string  =  "C:/Users/VERBISTF/OneDrive - KU Leuven/PhD Flore/P4_negative_emissions/Figures/Plots/boxplot_$(i).svg"
        # savefig(filename_string) 
        # savefig(p, filename_string)
        p = 0
    end



end

function scenario_bars_unfolded!(data)
    Sensitivity_tab = CSV.read("Sensitivity_dict_table.csv", DataFrame)
    Sensitivity_tab = data_original # This contains also the no removals data - as it was added later --> see processing function
    Sensitivity_tab[!, 2] = [eval(Meta.parse(replace(Sensitivity_tab[i,2], "Any" => ""))) for i in 1:length(Sensitivity_tab[!,1])] 
    data = Sensitivity_tab 
    Case_data = CSV.read("Table_for_paper_sept.csv", DataFrame)
    #Case_data[:, "Total_budget_per_net"] = Case_data[:, "Total_budget"]
    variable_names = ["Total_net_emissions", "Total_emissions", "Emitter_cost_avg", "Total_budget_per_net",  "Learning_gains_per_rem", "Total_costs_per_net"]
    plot_title = ["Total net emissions", "Total gross emissions", "Emitter balance", "CCB balance", "Learning gains removals", "Total cost per net emission"]
    plot_title = ["","","","","","","","",""]
    y_title = [L"Net emissions [GtCO$_2$]",  L"Gross emissions [GtCO$_2$]", L"Emitters' costs per gross emission [€$$/tCO$_2$]", L"Income per net emission [€$$/tCO$_2$]", L"Economic gains per gross removal [€$$/tCO$_2$]",  L"Total costs per net emission [€$$/tCO$_2$]"]
    Model_names = [L"\parbox{7em}{\centering Social \\ Optimum}", L"\parbox{7em}{\centering No \\ Removals}", L"\parbox{7em}{\centering Target-based}", L"\parbox{7em}{\centering Net Emission \\ Cap}",    L"\parbox{7em}{\centering Removal \\ Obligation}", L"\parbox{7em}{\centering Equivalence \\ Obligation}",  L"\parbox{7em}{\centering Conditional \\ Integration}",  L"\parbox{7em}{\centering Gross Emission \\ Cap}",   L"\parbox{7em}{\centering Gross cap \\ Obligation}"]
    MC = ["springgreen4", "black", "steelblue4", "deepskyblue3", "slateblue3", "thistle", "mediumorchid2", "chocolate1", "brown3"]
    ordered_Cnr = [2, 3, 4, 7, 8, 6, 5, 9]
    scattercolors = ["black" "gray73" "white"] #range(colorant"red", stop=colorant"green", length=length(data[1][1]))
    scatteralpha = [0.9 0.5 0.1]
    scattermarkershape = [:circle, :rect, :utriangle]
    scattersize = [2,3,4]
    markerwidth = [0.5, 1.5]
    ymin_vec = [10, 10, 100, -500, -20, 0]
    ymax_vec = [35, 60, 800, 800, 250, 700]
    # index_nr = 6
    # i = variable_names[index_nr]
    for (index_nr,i) in enumerate(variable_names)

        ymax = maximum(hcat([data[data[!,1] .== "(\"$(CASE_NAME[n])\", \"$(i)\")", 2][1][1:end] for n in ordered_Cnr ]...))
        ymin = minimum(hcat([data[data[!,1] .== "(\"$(CASE_NAME[n])\", \"$(i)\")", 2][1][1:end] for n in ordered_Cnr ]...))
        ydelta = 0.05 * (ymax - ymin)
        YMAX = ymax + ydelta
        YMIN = ymin - ydelta

    

        p = boxplot(reshape(hcat([data[data[!,1] .== "(\"$(CASE_NAME[n])\", \"$(i)\")", 2][1][1:end] for n in ordered_Cnr ]...), 54,length(ordered_Cnr)), xticks = ([m for m in 1:length(ordered_Cnr)], [Model_names[n] for n in 2:(1+length(ordered_Cnr))]), alpha = 0.5, label = "", title = "$(plot_title[index_nr])", ylabel = y_title[index_nr], xlabel = "", legend = false, xrotation=90, color = hcat([MC[m] for m in 2:(1+length(ordered_Cnr))]...), ylim = (YMIN, YMAX), xtickfontsize=12,ytickfontsize=12, xguidefontsize=14,yguidefontsize=14)
      
        if index_nr .== 1
            
            for i in 1:1
                annotate!([0.7], [19], text(L"1", 10))
                annotate!([0.85], [19], text(L"2", 10))
                annotate!([1.15], [19], text("3", 10))
                annotate!([1.3], [19], text("4", 10))
            end
        else 
            skip 
        end
        # for j in 1:1:length(data[!,2][1])
        #     scatter!([m for m in 1:length(ordered_Cnr)],[data[data[!,1] .== "(\"$(CASE_NAME[n])\", \"$(i)\")", 2][1][j] for n in ordered_Cnr ],  color = [scattercolors[M_r_vector[j]] for m in 1:length(ordered_Cnr)], markershape = [scattermarkershape[D_r_vector[j]] for m in 1:length(ordered_Cnr)], markersize = [scattersize[D_a_vector[j]] for m in 1:length(ordered_Cnr)], markerstrokewidth =[markerwidth[Equivalence_vector[j]] for m in 1:length(ordered_Cnr)], ylim = (YMIN, YMAX))
        # end
        Mr_2 = [x == 2 ? 1 : 0 for x in M_r_vector]
        Dr_2 = [x == 2 ? 1 : 0 for x in D_r_vector]
        Da_2 = [x == 2 ? 1 : 0 for x in D_a_vector]
        EQ_yes = [x == 1 ? 1 : 0 for x in Equivalence_vector]
        range = [i for i in 1:1:length(data[!,2][1])]
        range_Mr = [range[x] for x in 1:length(range) if (Dr_2[x] == 1) && (Da_2[x] == 1) && (EQ_yes[x] == 1)][[1,3]]
        range_Dr = [range[x] for x in 1:length(range) if (Mr_2[x] == 1) && (Da_2[x] == 1) && (EQ_yes[x] == 1)][[1,3]]
        range_Da = [range[x] for x in 1:length(range) if (Dr_2[x] == 1) && (Mr_2[x] == 1) && (EQ_yes[x] == 1)][[1,3]]
        range_EQ = [range[x] for x in 1:length(range) if (Dr_2[x] == 1) && (Da_2[x] == 1) && (Mr_2[x] == 1)][2]
        


        for j in range_Mr
            # effect of variations in removal learning (Mr) 
            scatter!([m-0.3 for m in 1:length(ordered_Cnr)],[data[data[!,1] .== "(\"$(CASE_NAME[n])\", \"$(i)\")", 2][1][j] for n in ordered_Cnr ],  color = [scattercolors[M_r_vector[j]] for m in 1:length(ordered_Cnr)], markershape = [scattermarkershape[D_r_vector[j]] for m in 1:length(ordered_Cnr)], markersize = [scattersize[D_a_vector[j]] for m in 1:length(ordered_Cnr)], markerstrokewidth =[markerwidth[Equivalence_vector[j]] for m in 1:length(ordered_Cnr)], ylim = (YMIN, YMAX))
        end

        for j in range_Dr
            # effect of variations in removal cost (Dr) 
            scatter!([m-0.15 for m in 1:length(ordered_Cnr)],[data[data[!,1] .== "(\"$(CASE_NAME[n])\", \"$(i)\")", 2][1][j] for n in ordered_Cnr ],  color = [scattercolors[M_r_vector[j]] for m in 1:length(ordered_Cnr)], markershape = [scattermarkershape[D_r_vector[j]] for m in 1:length(ordered_Cnr)], markersize = [scattersize[D_a_vector[j]] for m in 1:length(ordered_Cnr)], markerstrokewidth =[markerwidth[Equivalence_vector[j]] for m in 1:length(ordered_Cnr)], ylim = (YMIN, YMAX))
        end

        for j in range_Da
            # effect of variations in abatement cost (Da) 
            scatter!([m+0.15 for m in 1:length(ordered_Cnr)],[data[data[!,1] .== "(\"$(CASE_NAME[n])\", \"$(i)\")", 2][1][j] for n in ordered_Cnr ],  color = [scattercolors[M_r_vector[j]] for m in 1:length(ordered_Cnr)], markershape = [scattermarkershape[D_r_vector[j]] for m in 1:length(ordered_Cnr)], markersize = [scattersize[D_a_vector[j]] for m in 1:length(ordered_Cnr)], markerstrokewidth =[markerwidth[Equivalence_vector[j]] for m in 1:length(ordered_Cnr)], ylim = (YMIN, YMAX))

        end

        for j in range_EQ
            # effect of variations in equivalence factor (phi) 
            scatter!([m+0.3 for m in 1:length(ordered_Cnr)],[data[data[!,1] .== "(\"$(CASE_NAME[n])\", \"$(i)\")", 2][1][j] for n in ordered_Cnr ],  color = [scattercolors[M_r_vector[j]] for m in 1:length(ordered_Cnr)], markershape = [scattermarkershape[D_r_vector[j]] for m in 1:length(ordered_Cnr)], markersize = [scattersize[D_a_vector[j]] for m in 1:length(ordered_Cnr)], markerstrokewidth =[markerwidth[Equivalence_vector[j]] for m in 1:length(ordered_Cnr)], ylim = (YMIN, YMAX))
        end

        for j in [27]
            # effect of variations in removal learning (Mr) 
            scatter!([m for m in 1:length(ordered_Cnr)],[data[data[!,1] .== "(\"$(CASE_NAME[n])\", \"$(i)\")", 2][1][j] for n in ordered_Cnr ],  color = [scattercolors[M_r_vector[j]] for m in 1:length(ordered_Cnr)], markershape = [scattermarkershape[D_r_vector[j]] for m in 1:length(ordered_Cnr)], markersize = [scattersize[D_a_vector[j]] for m in 1:length(ordered_Cnr)], markerstrokewidth =[markerwidth[Equivalence_vector[j]] for m in 1:length(ordered_Cnr)], ylim = (YMIN, YMAX))
        end

        # here note that the actual base run generates the same results as sensitivity analysis, but for the equivalence case, the base run calculates with a 0.75 permanence fraction, while in the sensitivity runs, the impermancence is determined by the permanence in advance, using the same fraction
        # This means that the marker will be higher than the medium-grey square.
        scatter!([m for m in 1:length(ordered_Cnr)], [Case_data[(m+1), i] for m in 1:length(ordered_Cnr)], markershape = [:xcross], markeralpha = [0.5], markersize = [8], markercolor = [:black], markerstrokewidth =[3])
        # for j in [27]
        #     # effect of variations in removal learning (Mr) 
        #     scatter!([m for m in 1:length(ordered_Cnr)],[data[data[!,1] .== "(\"$(CASE_NAME[n])\", \"$(i)\")", 2][1][j] for n in ordered_Cnr ],  color = [scattercolors[M_r_vector[j]] for m in 1:length(ordered_Cnr)], markershape = [:xcross], markeralpha = [0.5], markersize = [8], markercolor = [:black], markerstrokewidth =[3])
        # end


        plot!(twinx(p), label="",  color = :white, xformatter=_->"", yformatter=_->"", ylim = (YMIN, YMAX)) 
        plot!(twiny(p), label="",  color = :white, xformatter=_->"", yformatter=_->"", gridlinewidth = 0, grid =:hide, gridalpha = 0, gridstyle = :dash, ylim = (YMIN, YMAX)) 
        plot!([2.5, 2.5], [ymin_vec[index_nr], ymax_vec[index_nr]], ls = :dash, color = :black, label = "")
        plot!([6.5, 6.5], [ymin_vec[index_nr], ymax_vec[index_nr]], ls = :dash, color = :black, label = "")
        plot!(grid = true, size  = (600, 400), ylim = (ymin_vec[index_nr], ymax_vec[index_nr])) # ylim = (0, 700)
        if index_nr == 1


            k = 0.2
            q = 0

            # Removal learning (Mr)
            annotate!(2.2+k, 34+q, text(L"\parbox{20em}{1: Removal learning: \quad   (low), \quad  (high)}", 8,:white, :black))
            scatter!([2.25+k],[34+q],  color = [scattercolors[M_r_vector[range_Mr[1]]] for m in 1:length(ordered_Cnr)], markershape = [scattermarkershape[D_r_vector[range_Mr[1]]] for m in 1:length(ordered_Cnr)], markersize = [scattersize[D_a_vector[range_Mr[1]]] for m in 1:length(ordered_Cnr)], markerstrokewidth =[markerwidth[Equivalence_vector[range_Mr[1]]] for m in 1:length(ordered_Cnr)], ylim = (YMIN, YMAX))
            scatter!([2.99+k],[34+q],  color = [scattercolors[M_r_vector[range_Mr[2]]] for m in 1:length(ordered_Cnr)], markershape = [scattermarkershape[D_r_vector[range_Mr[2]]] for m in 1:length(ordered_Cnr)], markersize = [scattersize[D_a_vector[range_Mr[2]]] for m in 1:length(ordered_Cnr)], markerstrokewidth =[markerwidth[Equivalence_vector[range_Mr[1]]] for m in 1:length(ordered_Cnr)], ylim = (YMIN, YMAX))

            # Removal cost (Dr)
            annotate!(2.2+k, 33+q, text(L"\parbox{20em}{2: Removal costs:\quad  \quad \quad (low),  \quad  (high)}", 8,:white, :black))
            scatter!([2.25+k],[33+q],  color = [scattercolors[M_r_vector[range_Dr[1]]] for m in 1:length(ordered_Cnr)], markershape = [scattermarkershape[D_r_vector[range_Dr[1]]] for m in 1:length(ordered_Cnr)], markersize = [scattersize[D_a_vector[range_Dr[1]]] for m in 1:length(ordered_Cnr)], markerstrokewidth =[markerwidth[Equivalence_vector[range_Dr[1]]] for m in 1:length(ordered_Cnr)], ylim = (YMIN, YMAX))
            scatter!([2.99+k],[33+q],  color = [scattercolors[M_r_vector[range_Dr[2]]] for m in 1:length(ordered_Cnr)], markershape = [scattermarkershape[D_r_vector[range_Dr[2]]] for m in 1:length(ordered_Cnr)], markersize = [scattersize[D_a_vector[range_Dr[2]]] for m in 1:length(ordered_Cnr)], markerstrokewidth =[markerwidth[Equivalence_vector[range_Dr[1]]] for m in 1:length(ordered_Cnr)], ylim = (YMIN, YMAX))


            # Abatement cost (Da)
            annotate!(2.2+k, 32+q, text(L"\parbox{20em}{3: Abatement costs:   \quad  \hspace{0.1mm}  (low), \quad   (high)}", 8,:white, :black))
            scatter!([2.25+k],[32+q],  color = [scattercolors[M_r_vector[range_Da[1]]] for m in 1:length(ordered_Cnr)], markershape = [scattermarkershape[D_r_vector[range_Da[1]]] for m in 1:length(ordered_Cnr)], markersize = [scattersize[D_a_vector[range_Da[1]]] for m in 1:length(ordered_Cnr)], markerstrokewidth =[markerwidth[Equivalence_vector[range_Da[1]]] for m in 1:length(ordered_Cnr)], ylim = (YMIN, YMAX))
            scatter!([2.99+k],[32+q],  color = [scattercolors[M_r_vector[range_Da[2]]] for m in 1:length(ordered_Cnr)], markershape = [scattermarkershape[D_r_vector[range_Da[2]]] for m in 1:length(ordered_Cnr)], markersize = [scattersize[D_a_vector[range_Da[2]]] for m in 1:length(ordered_Cnr)], markerstrokewidth =[markerwidth[Equivalence_vector[range_Da[1]]] for m in 1:length(ordered_Cnr)], ylim = (YMIN, YMAX))


            # Permanence factor (phi)
            annotate!(2.2+k, 31+q, text(L"\parbox{20em}{4: Permanence: \quad \quad \quad  \hspace{0.5mm} (no)}", 8,:white, :black))
            scatter!([2.25+k],[31+q],  color = [scattercolors[M_r_vector[range_EQ[1]]] for m in 1:length(ordered_Cnr)], markershape = [scattermarkershape[D_r_vector[range_EQ[1]]] for m in 1:length(ordered_Cnr)], markersize = [scattersize[D_a_vector[range_EQ[1]]] for m in 1:length(ordered_Cnr)], markerstrokewidth =[markerwidth[Equivalence_vector[range_EQ[1]]] for m in 1:length(ordered_Cnr)], ylim = (YMIN, YMAX))

            # Baserun 
            annotate!(2.2+k, 30+q, text(L"\parbox{20em}{ \quad \hspace{0.1mm} Base run: \quad \quad \quad  \quad \hspace{0.3mm}  / }", 8,:white, :black))
            scatter!([2.55+k], [30+q], markershape = [:xcross], markeralpha = [0.5], markersize = [8], markercolor = [:black], markerstrokewidth =[3])
            scatter!([2.25+k],[30+q],  color = [scattercolors[M_r_vector[27]] for m in 1:length(ordered_Cnr)], markershape = [scattermarkershape[D_r_vector[27]] for m in 1:length(ordered_Cnr)], markersize = [scattersize[D_a_vector[27]] for m in 1:length(ordered_Cnr)], markerstrokewidth =[markerwidth[Equivalence_vector[27]] for m in 1:length(ordered_Cnr)], ylim = (YMIN, YMAX))

        else 
            skip 
        end
        
        display(p)
        filename_string  =  "C:/Users/VERBISTF/OneDrive - KU Leuven/PhD Flore/P4_negative_emissions/Figures/Plots/boxplot_$(i).svg"
        savefig(filename_string) 
        savefig(p, filename_string)
        p = 0
    end



end


function ETS_variance!(ETS_data::Dict)
Model_names = [L"\parbox{20em}{\centering 1.1 NEM cap}", L"\parbox{20em}{\centering 1.2 NEM cap conditional}", L"\parbox{20em}{\centering 2.1 NAL cap}", L"\parbox{20em}{\centering 2.2 NAL \& NEM}"]
Key_names = ["NEM_cap" "NEM_conditional" "NAL_AdjCap" "NAL_cap"]

for i in 1:length(Model_names)
    Key = Key_names[i]
    Mean_data = [mean(ETS_data[(Key,t)]) for t in 1:T_OPT]
    Median_data =[median(ETS_data[(Key,t)]) for t in 1:T_OPT]
    p25 = [percentile(ETS_data[(Key,t)], 25) for t in 1:T_OPT]
    p75 = [percentile(ETS_data[(Key,t)], 75) for t in 1:T_OPT]
    max = [maximum(ETS_data[(Key,t)]) for t in 1:T_OPT]
    min = [minimum(ETS_data[(Key,t)]) for t in 1:T_OPT]
    x = Timeseries_data[1:1:(T_OPT),"Year"]

    plot(x, Mean_data, label="Mean", color = "black", linewidth =2.5)
    plot!(x, Median_data, label="Median", color = "grey", linewidth =2.5)
    plot!(x, p25, label="", color = "black", linestyle = :dash)
    plot!(x, p75, label= "",  color = "black", linestyle = :dash)
    plot!(x, min, label= "",  color = "black", linestyle = :dot)
    plot!(x, max, label= "Min-Max",  color = "black", linestyle = :dot)

    p= plot!(x, p25, fillrange = p75, color = "grey", fillalpha = 0.35, label = "Confidence band (p25-p75)", legend = :topleft, title = L"%$(Model_names[i])", ylabel = L"€/tCO$_2$", xlabel=L"Years", ylim=[0,750])
    display(p)
    filename_string  =  "C:/Users/VERBISTF/OneDrive - KU Leuven/PhD Flore/P4_negative_emissions/Figures/Plots/Variance_$(Key_names[i]).svg"
    savefig(filename_string) 
    savefig(p, filename_string)
    p = 0
end

end


function dict_table(Dict_results::Dict)


    
    
    Dict_table = DataFrame(Case_name = Plot_title_ordered, 
    Total_net_emissions = [Dict_results[i, "Total_net_emissions"][end] for i in CASE_NAME_ordered], 
    Total_emissions = [Dict_results[i, "Total_emissions"][end] for i in CASE_NAME_ordered], 
    Total_costs = [Dict_results[i, "Total_costs_per_net"][end].*Dict_results[i, "Total_net_emissions"][end] for i in CASE_NAME_ordered],
    Fraction_costs_abatement = [Dict_results[i, "Fraction_costs_abat"][end].*100 for i in CASE_NAME_ordered],
    Fraction_costs_removals = [Dict_results[i, "Fraction_costs_rem"][end].*100 for i in CASE_NAME_ordered],
    Fraction_removals = [(1.0-Dict_results[i, "Fraction_removals"][end])*100 for i in CASE_NAME_ordered], 
    Total_costs_per_net = [Dict_results[i, "Total_costs_per_net"][end] for i in CASE_NAME_ordered],
    Total_budget_per_net = [Dict_results[i, "Total_budget_per_net"][end] for i in CASE_NAME_ordered],
    Learning_gains_per_em = [Dict_results[i, "Learning_gains_per_em"][end] for i in CASE_NAME_ordered],
    Learning_gains_per_rem = [Dict_results[i, "Learning_gains_per_rem"][end] for i in CASE_NAME_ordered],
    ETS_max = [Dict_results[i, "ETS_max"][end] for i in CASE_NAME_ordered], 
    RC_max = [Dict_results[i, "RC_max"][end] for i in CASE_NAME_ordered],
    ETS_avg = [Dict_results[i, "ETS_avg"][end] for i in CASE_NAME_ordered], 
    Emitter_cost_avg = [Dict_results[i, "Emitter_cost_avg"][end] for i in CASE_NAME_ordered])

    CSV.write("Table_for_paper_sept.csv", Dict_table)
    return Dict_table

end


function learningcurves(results, ADMM)

    # first run a model to get the results, ADMM dataframes
    p = plot(title = "", gridlinewidth = 1, xlims= (0, 1000), ylim=(0, 1250))
    shades_of_green = range(colorant"darkseagreen1", stop=colorant"green", length=length(1:T_OPT))
    st = 0.001
    for t in 2:T_OPT #model_r.ext[:sets][:T]
        st = 0.001
        rem_max = U_R./(1+exp(K_R*(t)))
        remov = collect(0:st:rem_max)
        rem_costs = [minimum([1000, fr(remov, t, results["Cum_remov"][ADMM["n_iter"],:])[i]]) for i in 1:length(remov)]
        plot!(remov*EM_REF, rem_costs, color = shades_of_green[t],  ylabel = L"Marginal Removal/Abatement Costs [€/tCO$_2$]", xlabel = L"Amount of removals/emissions [MtCO$_2$]", label = "")
    end
    rem_max = U_R./(1+exp(K_R*(T_OPT)))
    remov = collect(0:st:rem_max)
    rem_costs = [minimum([1000, fr(remov, T_OPT, results["Cum_remov"][ADMM["n_iter"],:])[i]]) for i in 1:length(remov)]
    plot!(remov*EM_REF, rem_costs, color = shades_of_green[T_OPT], label = "MRCC")

    shades_of_blue = range(colorant"steelblue", stop=colorant"navyblue", length=length(1:T_OPT))
    #Cum_abat_graph = Dict(i => results["Cum_abat"][ADMM["n_iter"],i+1] for i in 0:T_OPT )
    for t in 2:T_OPT #model_r.ext[:sets][:T]
        st = 0.001
        abat_max = U_A./(1+exp(K_A*(t)))
        abat = collect(0:st:abat_max)
        abat_costs = fa(abat, t, results["Cum_abat"][ADMM["n_iter"],:])
        plot!(abat*EM_REF, reverse(abat_costs), color = shades_of_blue[t],  xlabel = L"Amount of removals/emissions [MtCO$_2$]", label = "")
    end
    abat_max = U_A./(1+exp(K_A*(T_OPT)))
    abat = collect(0:st:abat_max)
    abat_costs = fa(abat, T_OPT, results["Cum_abat"][ADMM["n_iter"],:])
    plot!(abat*EM_REF, reverse(abat_costs), color = shades_of_blue[1], label = "MACC", ylims=(0,1300), xlims = (0,1000))
    annotate!(35, 750, Plots.text(L"\textbf{t$_0$}", 18, color = "darkgreen"), label=false)
    annotate!(50, 50, Plots.text(L"\textbf{t$_{end}$}", 18, color = "darkgreen"), label=false)
    annotate!(965, 50, Plots.text(L"\textbf{t$_{end}$}", 18, color = "dodgerblue4"), label=false)
    annotate!(950, 220, Plots.text(L"\textbf{t$_0$}", 18, color = "dodgerblue4"), label=false)
    plot!(twinx(p), [0 1000], [0 0], xlims= (0, 1000), ylim=(0, 1250), xticks = 0:200:1000, yticks = 0:200:1250, label="",  color = :black, xformatter=_->"", yformatter=_->"") 
    plot!(twiny(p), [0 1000], [0 0], xlims= (0, 1000), ylim=(0, 1250),  xticks = 0:200:1000, yticks = 0:200:1250, label="", color = :black, xformatter=_->"", yformatter=_->"")
    plot!(grid = true, gridalpha=0.03, gridlinewidth=1, xlims= (0, 1000), ylim=(0, 800),  legend = (0.80,0.97), xticks = 0:200:1000, yticks = 0:200:1250, xtickfontsize=12,ytickfontsize=12, legendfontsize = 14, yguidefontsize = 14, xguidefontsize = 14)


    # em_max = 1
    # em = (0:0.001:em_max)*EM_REF
    # f_a = [BETA.*(EM_REF - em[i]).^GAMMA for i in 1:length(em)]
    # plot!(em, f_a, color = :deepskyblue2, label = "MACC (noLearning)", xlim=(0,1000))
    # rem = (0:0.001:em_max)*EM_REF
    # f_r = [A_r + B_r*rem[i] for i in 1:length(rem)]
    # plot!(em, f_a, color = :deepskyblue2, label = "MACC (noLearning)", xlim=(0,1000))




    filename_string  =  "C:/Users/VERBISTF/OneDrive - KU Leuven/PhD Flore/P4_negative_emissions/Figures/Plots/Learning_curves_r2.svg"
    savefig(filename_string) 
    savefig(p, filename_string)
return  display(p)
end


function learning_gains_case(C_nr_1, C_nr_2)

    ADMM_max_it = 500
    Case_name_1 = CASE_NAME[C_nr_1]
    RHO = RHO1 = RHO2 = RHO_vect[C_nr_1]
    ###--------------------------------------------------------------
    results_1 = Dict()
    ADMM_1 = Dict()
    global BUDGET = Carbon_budget

    print("$(Case_name_1)  --  ")
    include("$(Case_name_1).jl")

    define_results!(Parameters,results_1,ADMM_1)

    model_e_1 = Model(optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" =>  0))
    model_r_1 = Model(optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" =>  0))
    model_ccb_1 = Model(optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" =>  0))


    # Step 3b: process parameters
    process_sets_parameters!(model_e_1, Parameters,Timeseries_data, T_OPT)
    process_sets_parameters!(model_r_1, Parameters,Timeseries_data, T_OPT)
    process_sets_parameters!(model_ccb_1, Parameters,Timeseries_data, T_OPT)


    # calculate equilibrium 
    ADMM_1["walltime"] = @elapsed Dict_results = ADMM_price_based_m!(Case_name_1,results_1, ADMM_1, model_e_1, model_r_1, model_ccb_1, Dict_results, LINEAR)  

    ADMM_max_it = 500
    Case_name_2 = CASE_NAME[C_nr_2]
    RHO = RHO1 = RHO2 = RHO_vect[C_nr_2]
    ###--------------------------------------------------------------
    results_2 = Dict()
    ADMM_2 = Dict()
    global BUDGET = Carbon_budget

    print("$(Case_name_2)  --  ")
    include("$(Case_name_2).jl")

    define_results!(Parameters,results_2,ADMM_2)

    model_e_2 = Model(optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" =>  0))
    model_r_2 = Model(optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" =>  0))
    model_ccb_2 = Model(optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" =>  0))


    # Step 3b: process parameters
    process_sets_parameters!(model_e_2, Parameters,Timeseries_data, T_OPT)
    process_sets_parameters!(model_r_2, Parameters,Timeseries_data, T_OPT)
    process_sets_parameters!(model_ccb_2, Parameters,Timeseries_data, T_OPT)


    # calculate equilibrium 
    ADMM_2["walltime"] = @elapsed Dict_results = ADMM_price_based_m!(Case_name_2,results_2, ADMM_2, model_e_2, model_r_2, model_ccb_2, Dict_results, LINEAR)  



    parameters_extract(model_r_1)
    Cum_remov_extract_1 = Dict( i =>  results_1["Cum_remov"][ADMM_1["n_iter"],i+1] for i in 0:T_OPT)
    Cum_remov_extract_1[-1] = Cum_remov_init # creation of element in case that t = 0
    L_1 = [1/(1+DRc)^t*quadgk(var -> f(var, t, Cum_remov_extract_1), 0, r[t]./EM_REF)[1].*EM_REF for t in model_r_1.ext[:sets][:T]]
    NL_1 = [1/(1+DRc)^t*quadgk(var -> f(var, t, Cum_NoLearning_remov), 0, r[t]./EM_REF)[1].*EM_REF for t in model_r_1.ext[:sets][:T]]
    L_1_vector = [L_1[i][1] for i in 1:length(L_1)]
    NL_1_vector = [NL_1[i][1] for i in 1:length(L_1)]
    R1_vector = [r[i] for i in 0:(length(r)-1)]
    plot(L_1_vector./R1_vector, label = "Costs removal - learning $(Case_name_1)") 
    plot!(NL_1_vector./R1_vector, label = "Costs removal - no learning $(Case_name_1)") 


    parameters_extract(model_r_2)
    Cum_remov_extract_2 = Dict( i =>  results_2["Cum_remov"][ADMM_2["n_iter"],i+1] for i in 0:T_OPT)
    Cum_remov_extract_2[-1] = Cum_remov_init # creation of element in case that t = 0
    L_2 = [1/(1+DRc)^t*quadgk(var -> f(var, t, Cum_remov_extract_2), 0, r[t]./EM_REF)[1].*EM_REF for t in model_r_2.ext[:sets][:T]]
    NL_2 = [1/(1+DRc)^t*quadgk(var -> f(var, t, Cum_NoLearning_remov), 0, r[t]./EM_REF)[1].*EM_REF for t in model_r_2.ext[:sets][:T]]
    L_2_vector = [L_2[i][1] for i in 1:length(L_2)]
    NL_2_vector = [NL_2[i][1] for i in 1:length(L_2)]
    R2_vector = [r[i] for i in 0:(length(r)-1)]
    plot!(L_2_vector./R2_vector, label = "Costs removal - learning $(Case_name_2)") 
    plot!(NL_2_vector./R2_vector, label = "Costs removal - no learning $(Case_name_2)") 


end

function deflation_effect()

    ADMM_max_it = ADMM_max_it
    C_nr_1 = C_nr_2 =  8  # equivalence obligations
    Case_name_1 = CASE_NAME[C_nr_1]
    RHO = RHO1 = RHO2 = RHO_vect[C_nr_1]
    FRAC = -0.75
    ###--------------------------------------------------------------
    results_1 = Dict()
    ADMM_1 = Dict()
    global BUDGET = Carbon_budget

    print("$(Case_name_1)  --  ")
    include("$(Case_name_1).jl")

    define_results!(Parameters,results_1,ADMM_1)

    model_e_1 = Model(optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" =>  0))
    model_r_1 = Model(optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" =>  0))
    model_ccb_1 = Model(optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" =>  0))


    # Step 3b: process parameters
    process_sets_parameters!(model_e_1, Parameters,Timeseries_data, T_OPT)
    process_sets_parameters!(model_r_1, Parameters,Timeseries_data, T_OPT)
    process_sets_parameters!(model_ccb_1, Parameters,Timeseries_data, T_OPT)


    # calculate equilibrium 
    ADMM_1["walltime"] = @elapsed Dict_results = ADMM_price_based_m!(Case_name_1,results_1, ADMM_1, model_e_1, model_r_1, model_ccb_1, Dict_results, LINEAR)  

    #---------
    parameters_extract(model_e_1)
    parameters_extract(model_r_1)


    # Step 6: 
    EM_1 = results_1["e"][ADMM_1["n_iter"],:]
    REM_1 = results_1["r"][ADMM_1["n_iter"],:]
    hatch_vect = [x for p in zip(EM, zeros(length(EM))) for x in p].*[0; repeat([1, 1, 0, 0], convert(Int64, T_OPT/2)); 1]
    Doubled_time_vector =  [x for x in Timeseries_data[1:1:(T_OPT+1),"Year"] for _ in 1:2] + [repeat([0, 0, 0, -1], convert(Int64, T_OPT/2)); 0; 0]

    Peac_1 = results_1["EAC_price"][ADMM_1["n_iter"],:]
    Pcdr_1 = results_1["CDR_price"][ADMM_1["n_iter"],:]

    variables_extract(model_e_1)
    variables_extract(model_r_1)
    variables_extract(model_ccb_1)
    
    # parameters_extract(model)
    # variables_extract(model)


    EAC_surplus = [sum( ec[i]  - (EM[i+1])  for i in 0:t) for t in model_r_1.ext[:sets][:T]]
    CDR_surplus = [sum((REM[i+1] -rc[i])  for i in 0:t) for t in model_r_1.ext[:sets][:T]]
    CCB_surplus = [sum((rc_r[i] -rc_e[i])  for i in 0:t) for t in model_r_1.ext[:sets][:T]]



    ymax = 500
    Color = ["snow3" "snow2"  "slategray3"]
    Alpha = [1 0.4 0.5]

    ppp = groupedbar(Timeseries_data[1:1:(T_OPT+1),"Year"], [EM_1 -REM_1], bar_position =:stack, title = "Deflation Effect", xticks = (2025:5:(2024+T_OPT)),  labels = ["Emissions"  "Removals"], ylabel = L"Emissions/Removals [MtCO$_2$/yr]", xlabel = "Years", ylims = (-600, 2500), yticks =-500:250:2500, foreground_color_legend = nothing,  legend = (0.02, 0.88), lw = 0.75, color = Color, linecolor = [:gray31 :gray31], alpha = [1.0 0.8], right_margin = 5mm)
    plot!(twiny(ppp), Timeseries_data[1:1:(T_OPT+1),"Year"], [Timeseries_data[1:1:(T_OPT+1), "ETS_cap [EUA/y]"]./1], lw=[1 1], xticks = (2025:5:(2024+T_OPT)), ylims = (-600, 2500), label = ["ETS cap"],foreground_color_legend = nothing, xformatter=_->"", yformatter=_->"", xlabel = "", ylabel="", legend = (0.02,0.98), color = ["black"  "black" "black" "brown"], ls = [:solid :dash :dash :dash], margin = 0mm)
    plot!(twinx(ppp), Timeseries_data[1:1:(T_OPT+1),"Year"], [[Peac_1[t]   for t in 1:1:(T_OPT+1)] [Pcdr_1[t]  for t in 1:1:(T_OPT+1)]], ylims = (0, ymax), yticks = 0:50:ymax, foreground_color_legend = nothing, label = ["ETS price (MAC)" "RC price (MRC)"],  color = ["grey" "grey83"], legend = (0.65,0.98), lw = 2, margin = 0mm)  
    plot!(twiny(ppp), Timeseries_data[1:1:(T_OPT+1),"Year"], [(EM_1-REM_1)], lw=0, xticks = (2025:5:(2024+T_OPT)), ylims = (-600, 2500), label = ["Net emissions"],foreground_color_legend = nothing, xformatter=_->"", yformatter=_->"", xlabel = "", ylabel="", legend = (0.02,0.93), color = ["black"  "black" "black" "brown"], ls = [:dash :dash :dash :dash], margin = 0mm, markershape = :xcross, markersize = 2)

    #pprint = plot(ppp, size  = (500, 400),  xlim = [2023,(T_OPT_plot)],  titlefont = font(11),  title = Plot_title[C_nr], plot_titlefontsize = 12, topmargin = 2mm, gridlinewidth = 1)
    


    ###--------------------------------------------------------------
    
    #----------------------------------------
    ADMM_max_it = ADMM_max_it
    Case_name_2 = CASE_NAME[C_nr_2]
    RHO = RHO1 = RHO2 = RHO_vect[C_nr_2]
    FRAC = 0.0
    results_2 = Dict()
    ADMM_2 = Dict()
    global BUDGET = Carbon_budget

    print("$(Case_name_2)  --  ")
    include("$(Case_name_2).jl")

    define_results!(Parameters,results_2,ADMM_2)

    model_e_2 = Model(optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" =>  0))
    model_r_2 = Model(optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" =>  0))
    model_ccb_2 = Model(optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" =>  0))


    # Step 3b: process parameters
    process_sets_parameters!(model_e_2, Parameters,Timeseries_data, T_OPT)
    process_sets_parameters!(model_r_2, Parameters,Timeseries_data, T_OPT)
    process_sets_parameters!(model_ccb_2, Parameters,Timeseries_data, T_OPT)


    # calculate equilibrium 
    ADMM_2["walltime"] = @elapsed Dict_results = ADMM_price_based_m!(Case_name_2,results_2, ADMM_2, model_e_2, model_r_2, model_ccb_2, Dict_results, LINEAR)  


    parameters_extract(model_e_2)
    parameters_extract(model_r_2)


    # Step 6: 
    EM = results_2["e"][ADMM_2["n_iter"],:]
    REM = results_2["r"][ADMM_2["n_iter"],:]
    hatch_vect = [x for p in zip(EM, zeros(length(EM))) for x in p].*[0; repeat([1, 1, 0, 0], convert(Int64, T_OPT/2)); 1]
    Doubled_time_vector =  [x for x in Timeseries_data[1:1:(T_OPT+1),"Year"] for _ in 1:2] + [repeat([0, 0, 0, -1], convert(Int64, T_OPT/2)); 0; 0]

    Peac = results_2["EAC_price"][ADMM_2["n_iter"],:]
    Pcdr = results_2["CDR_price"][ADMM_2["n_iter"],:]

    variables_extract(model_e_2)
    variables_extract(model_r_2)
    variables_extract(model_ccb_2)
    
    # parameters_extract(model)
    # variables_extract(model)


    EAC_surplus = [sum( ec[i]  - (EM[i+1])  for i in 0:t) for t in model_r_2.ext[:sets][:T]]
    CDR_surplus = [sum((REM[i+1] -rc[i])  for i in 0:t) for t in model_r_2.ext[:sets][:T]]
    CCB_surplus = [sum((rc_r[i] -rc_e[i])  for i in 0:t) for t in model_r_2.ext[:sets][:T]]



    ymax = 500
    Color = ["snow3" "snow2"  "slategray3"]
    Alpha = [1 0.4 0.5]

    groupedbar!(ppp, Timeseries_data[1:1:(T_OPT+1),"Year"], [EM -REM], bar_position =:stack, label = nothing, xticks = (2025:5:(2024+T_OPT)), foreground_color_legend = nothing,   lw = 0.75, color = Color, linecolor = [:gray31 :gray31], ls = [:dash :dash], alpha = Alpha, right_margin = 5mm)
    groupedbar!(ppp, Timeseries_data[1:1:(T_OPT+1),"Year"], [EM_1 -REM_1], bar_position =:stack, xticks = (2025:5:(2024+T_OPT)),  label = nothing,  foreground_color_legend = nothing,  lw = 0.75, color = Color, linecolor = [:gray31 :gray31], alpha = [0.0 0.8], right_margin = 5mm)
    plot!(twinx(ppp), Timeseries_data[1:1:(T_OPT+1),"Year"], [[Peac[t]   for t in 1:1:(T_OPT+1)] [Pcdr[t]  for t in 1:1:(T_OPT+1)]], label = nothing, legend = false, alpha = [0.5 0.5],  ylims = (0, ymax), yticks = 0:50:ymax, foreground_color_legend = nothing, color = ["grey" "grey83"], ylabel = [L"Certificate price [€/tCO$_2$]"], lw = 2, margin = 0mm, ls = [:dash :dash])  

    pprint = plot(ppp, size  = (500, 400),  xlim = [2023,(T_OPT_plot)],  titlefont = font(11),  plot_titlefontsize = 12, topmargin = 2mm, gridlinewidth = 1)
    
    FRAC = 0.15
    display(pprint)


    filename_string  =  "C:/Users/VERBISTF/OneDrive - KU Leuven/PhD Flore/P4_negative_emissions/Figures/Plots/deflation_effect.svg"
    savefig(filename_string) 
    savefig(pprint, filename_string)

end


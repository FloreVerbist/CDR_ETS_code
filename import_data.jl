function import_data_CCUS(file::Any)
    #file ="C:/Users/VERBISTF/OneDrive - KU Leuven/Code examples/ets-ncc-master/BECCS/Script_BECCS/BECCS_data.xlsx"
    data = XLSX.readxlsx(file)
    Parameters = DataFrame(XLSX.readtable(file, "Parameters"))
    Parameters[!, :"Value"] = convert.(Float64, Parameters[!, :"Value"])
    Index_names = Parameters[!, "Parameter_name"]
    Timeseries_data = DataFrame(XLSX.readtable(file, "Timeseries"))
    return Timeseries_data, Index_names, Parameters
end


#Dictionaries for figure labels

function ranges_dict()
    d = Dict(
        "α"  => (-0.2, 1.2),
        "w0" => (-10, 90),
        "τ"  => (-0.2, 4.5),
        "PR" => (-0.2, 2.0),
        "u_ms"=>(-1, 11),
        "L_km"=>(100, 2100)
    )
    return d
end

function short_labels_dict()
    d = Dict{String, LaTeXString}(
        "spwp" => L"s_\mathrm{pwp}",
        "sfc"  => L"s_\mathrm{fc}",
        "ep"   => L"e_\mathrm{p}", 
        "eo"   => L"e_\mathrm{o}",
        "ϵ"    => L"\epsilon",
        "r"    => L"r",  
        "α"    => L"\alpha", 
        "nZr"  => L"nz_\mathrm{r}",   
        "a"    => L"a",  
        "b"    => L"b",
        "wsat" => L"w_\mathrm{sat}", 
        "u"    => L"u",
        "L"    => L"L",
        "Pl"   => L"P_\ell",
        "Po"   => L"$P_\mathrm{o}$",
        "El"   => L"$E_\ell$",
        "R"    => L"$R$",
        "PR"   => L"$\chi$",
        "Ptot" => L"$P_\mathrm{mean}$",
        "A"    => L"$A_\ell$",
        "B"    => L"$-A_\mathrm{o}$", 
        "τ"    => L"\tau",
        "τ_time"=> L"\tau^{-1}",
        "s"    => L"$s$",
        "wl"   => L"$w_\ell$",
        "wo"   => L"$w_\mathrm{o}$]",
        "L1_km"   => L"$L_1$",
        "L2_km"   => L"$L_2$",
        "L3_km"   => L"$L_3$",
        "u_ms"    => L"$u$",
        "w0"   => L"$w_0$",
        "PR12" => L"$PR_{12}$",
        "P1"   => L"$P_1$",
        "P2"   => L"$P_2$",
        "P3"   => L"$P_3$",
        "Δwtot"=> L"$Δw_\mathrm{tot}$",
    )
    return d
end

function labels_dict()
    d = Dict{String, LaTeXString}(
        "spwp" => L"$s_\mathrm{pwp}$",
        "ep"   => L"$e_\mathrm{p}$ [mm/day]", 
        "eo"   => L"$e_\mathrm{o}$ [mm/day]",
        "ϵ"    => L"$\epsilon$",
        "r"    => L"$r$",  
        "α"    => L"$\alpha$", 
        "nZr"  => L"$nZ_\mathrm{r}$ [mm]",   
        "a"    => L"$a$",  
        "b"    => L"$b$",
        "wsat" => L"$w_\mathrm{sat}$ [mm]", 
        "u"    => L"$u$ [m/s]",
        "L"    => L"$L$ [km]",
        "Pl"   => L"$P_\ell$ [mm/day]",
        "Po"   => L"$P_\mathrm{o}$ [mm/day]",
        "El"   => L"$E_\ell$ [mm/day]",
        "R"    => L"$R$ [mm/day]",
        "χ"    => L"$\chi$",
        "Φ"    => L"Infiltration $\Phi$",
        "Ptot" => L"$P_\mathrm{mean}$ [mm/day]",
        "A"    => L"$A_\ell$ [mm/day]",
        "B"    => L"$-A_\mathrm{o}$ [mm/day]", 
        "τ"    => L"$\tau$ [1/day]",
        "τ_time" => L"$\tau^{-1}$ [day]",
        "s"    => L"$s$",
        "wl"   => L"$w_\ell$ [mm]",
        "wo"   => L"$w_\mathrm{o}$ [mm]",
        "L1_km"   => L"$L_1$ [km]",
        "L2_km"   => L"$L_2$ [km]",
        "L3_km"   => L"$L_3$ [km]",
        "u_ms"    => L"$u$ [m/s]",
        "w0"   => L"$w_0$ [mm]",
        "PR12" => L"$PR_{12}$",
        "P1"   => L"$P_1$ [mm/day]",
        "P2"   => L"$P_2$ [mm/day]",
        "P3"   => L"$P_3$ [mm/day]",
        "Δwtot"=> L"$Δw_\mathrm{tot} [mm]$",
        "t"    => L"$t$ [day]",
    )
    return d
end

function labels_norm_dict()
    d = Dict{String, LaTeXString}(
        "spwp" => L"$s_\mathrm{pwp}$",
        "ep"   => L"$e_\mathrm{p}$", 
        "eo"   => L"$e_\mathrm{o}$",
        "ϵ"    => L"$\epsilon$",
        "r"    => L"$r$",  
        "α"    => L"$\alpha$", 
        "nZr"  => L"$nZ_\mathrm{r}$",   
        "a"    => L"$a$",  
        "b"    => L"$b$",
        "wsat" => L"$w_\mathrm{sat}$", 
        "u"    => L"$u$ [m/s]",
        "L"    => L"$L$ [km]",
        "Pl"   => L"$P_\ell$",
        "Po"   => L"$P_\mathrm{o}$",
        "El"   => L"$E_\ell$",
        "R"    => L"$R$",
        "PR"   => L"$\chi$",
        "Ptot" => L"$P_\mathrm{mean}$",
        "A"    => L"$A_\ell$",
        "B"    => L"$-A_\mathrm{o}$", 
        "τ"    => L"$\tau$",
        "τ_time" => L"$\tau^{-1}$",
        "s"    => L"$s$",
        "wl"   => L"$w_\ell$",
        "wo"   => L"$w_\mathrm{o}$",
        "L1_km"   => L"$L_1$",
        "L2_km"   => L"$L_2$",
        "L3_km"   => L"$L_3$",
        "u_ms"    => L"$u$",
        "w0"   => L"$w_0$",
        "PR12" => L"$PR_{12}$",
        "P1"   => L"$P_1$",
        "P2"   => L"$P_2$",
        "P3"   => L"$P_3$",
        "Δwtot"=> L"$Δw_\mathrm{tot}$",
    )
    return d
end


function full_labels_dict()
    d = Dict{String, LaTeXString}(
        "spwp" => L"Permanent wilting point $s_\mathrm{pwp}$",
        "α"    => L"Land fraction $\alpha$", 
        "r"    => L"Runoff exponent $r$", 
        "ϵ"    => L"Runoff parameter $\epsilon$",
        "a"    => L"Precipitation parameter $a$",
        "b"    => L"Precipitation parameter $b$",
        "wsat" => L"Saturation water vapor path $w_\mathrm{sat}$ [mm]",
        "eo"   => L"Ocean evaporation rate $E_\mathrm{o}$ [mm/day]",
        "Eo"   => L"Ocean evaporation $E_\mathrm{o}$ [mm/day]",
        "PR"   => L"Precipitation ratio $\chi$",
        "PRmean"=>L"Mean precipitation ratio $\chi_\mathrm{mean}$",
        "τ"    => L"Atmospheric transport parameter $\tau$ [1/day]",
        "u_ms" => L"Wind speed $u$ [m/s]",
        "τ_time" => L"Timescale of atmospheric transport $\tau$ [day]",
        "Atot" => L"Advection rate $A\, \,  \left[10^8\,\mathrm{mm}^2\mathrm{/day}\right]$",
        "s"    => L"Soil moisture saturation $s$",
        "El"   => L"Evapotranspiration $E_\ell$ [mm/day]",
        "L_km" => L"Full domain length $L$ [km]",
        "L1"   => L"First ocean length $L_1$ [km]",
        "L2"   => L"Land length $L_2$ [km]",
        "L3"   => L"Second ocean length $L_3$ [km]",
        "w0"   => L"Boundary water vapor path $w_0$ [mm]",
        "PR12" => L"Precipitation ratio $P_1/P_2$",
        "Δwtot"=> L"Total advection $Δw_\mathrm{tot}$ [mm]",
        "dw"   => L"Moisture difference $w_\mathrm{o}-w_\ell$ [mm]",
        "w1"   => L"Mean water vapor path $w_{\mathrm{o},1}$ [mm]",
        "w2"   => L"Mean water vapor path $w_{\ell}$ [mm]",
        "w3"   => L"Mean water vapor path $w_{\mathrm{o},2}$ [mm]",
        "mean_wo"=> L"Mean ocean water vapor path $w_\mathrm{o}$",
        "P1"   => L"Precipitation rate $P_1$ [mm/day]",
        "P2"   => L"Precipitation rate $P_2$ [mm/day]",
        "P3"   => L"Precipitation rate $P_3$ [mm/day]",
        "Po"   => L"Ocean precipitation $P_\mathrm{o}$ [mm/day]",
        "Pl"   => L"Land precipitation $P_\ell$ [mm/day]",
        "R"    => L"Runoff $R$ [mm/day]",
        "Φ"    => L"Infiltration function $\Phi$",
        "wl"   => L"Mean land water vapor path $w_\ell$ [mm]",
        "wo"   => L"Mean ocean water vapor path $w_\mathrm{o}$ [mm]",
        "ep"   => L"Potential ET $e_\mathrm{p}$ [mm/day]",
        "shat" => L"$\hat{s}$ indepdent of $s_\mathrm{pwp}$ and $s_\mathrm{fc}$",
        "Ehat" => L"Evapotranspiration $E$ [mm/day]",
        "A"    => L"Land advection rate $A$ [mm/day]",
        "Al"   => L"Land advection $A_\ell$ [mm/day]",
        "B"    => L"Ocean advection rate $B$ [mm/day]",
        "Ao"   => L"Ocean advection $A_\mathrm{o}$ [mm/day]",
        "Tsrf" => L"Surface temperature $T_\mathrm{srf}$ [K]",
        "t"    => L"Time $t$ [day]",
        "umax" => L"Maximum wind speed $u$ [m/s]",
        "watbal"=> L"$E_\ell + E_\mathrm{o} - P_\ell - P_\mathrm{o}$ [mm/day]",
        "wind" => L"Wind verlocity $u$ [m/s]",
    )
    return d
end


function titles_dict()
    d = Dict{String, String}(
        "spwp" => "Permanent wilting point",
        "ep"   => "Potential evapotranspiration", 
        "eo"   => "Mean ocean evaporation rate",
        "ϵ"    => "Runoff parameter",
        "r"    => "Runoff parameter",  
        "α"    => "Land fraction", 
        "nZr"  => "Reservoir depth",   
        "a"    => "Precipitation parameter",  
        "b"    => "Precipitation parameter",
        "wsat" => "Saturation water vapor path", 
        "u"    => "Wind speed",
        "L"    => "Full domain size",
        "Pl"   => "Land precipitation rate",
        "Po"   => "Ocean precipitation rate",
        "El"   => "Land evapotranspiration rate",
        "R"    => "Runoff rate",
        "PR"   => "Precipitation ratio",
        "Ptot" => "Mean precipitation rate",
        "A"    => "Land advection rate",
        "B"    => "Ocean advection rate", 
        "Φ"  => "Infiltration function",
        "s"    => "Rel. soil moisture saturation",
        "wl"   => "Mean land water vapour path",
        "wo"   => "Mean ocean water vapour path",
        "w0"   => "Boundary water vapor path",
        "L1_km"=> "First ocean length",
        "L2_km"=> "Land length",
        "L3_km"=> "Second ocean length",
        "PR12" => "PR without second ocean",
        "P1"   => "1st ocean precipitation rate",
        "P2"   => "Land precipitation rate",
        "P3"   => "2nd ocean precipitation rate",
        "Δwtot"=> "Total advection",
    )
    # d = Dict{String, LaTeXString}(
    #     "spwp" => LaTeXString("Permanent wilting point"),
    #     "ep"   => LaTeXString("Potential evapotranspiration"), 
    #     "eo"   => LaTeXString("Mean ocean evaporation rate"),
    #     "ϵ"    => LaTeXString("Runoff parameter"),
    #     "r"    => LaTeXString("Runoff parameter"),  
    #     "α"    => LaTeXString("Land fraction"), 
    #     "nZr"  => LaTeXString("Reservoir depth"),   
    #     "a"    => LaTeXString("Precipitation parameter"),  
    #     "b"    => LaTeXString("Precipitation parameter"),
    #     "wsat" => LaTeXString("Saturation water vapor path"), 
    #     "u"    => LaTeXString("Wind speed"),
    #     "L"    => LaTeXString("Full domain size"),
    #     "Pl"   => LaTeXString("Land precipitation rate"),
    #     "Po"   => LaTeXString("Ocean precipitation rate"),
    #     "El"   => LaTeXString("Land evapotranspiration rate"),
    #     "R"    => LaTeXString("Runoff rate"),
    #     "PR"   => LaTeXString("Precipitation ratio"),
    #     "Ptot" => LaTeXString("Mean precipitation rate"),
    #     "A"    => LaTeXString("Land advection rate"),
    #     "B"    => LaTeXString("Ocean advection rate"), 
    # )
    return d
end


function labels_lin_dict()
    d = Dict{String, LaTeXString}(
        "eo"   => L"Ocean evaporation $e_\mathrm{o}$ [mm/day]",
        "e"    => L"Land evaporation slope $e$ [mm/day]",
        "r"    => L"Runoff fraction slope $r$",  
        "α"    => L"Land fraction $\alpha$", 
        "p"    => L"Precipitation slope $p$ [1/day]",
        "Pl"   => L"Land precipitation $P_\ell$ [mm/day]",
        "Po"   => L"Ocean precipitation $P_\mathrm{o}$ [mm/day]",
        "El"   => L"Land evapotranspiration $E_\ell$ [mm/day]",
        "R"    => L"Runoff $R$ [mm/day]",
        "χ"   =>  L"Precipitation ratio $\chi$",
        "Φ"    => L"Infiltration $\Phi$",
        "Al"    => L"Land advection $A_\ell$ [mm/day]",
        "Ao"    => L"Ocean advection $-A_\mathrm{o}$ [mm/day]", 
        "τ"    => L"Transport parameter $\tau$ [1/day]",
        "τ_time" => L"$\tau^{-1}$ [day]",
        "s"    => L"Soil moisture saturation $s$",
        "wl"   => L"Land water vapor path $w_\ell$ [mm]",
        "wo"   => L"Ocean water vapor path $w_\mathrm{o}$ [mm]",
    )
    return d
end
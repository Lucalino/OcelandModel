# Functions used in the closed model analysis
# File created on 8th October 2021 by Luca Schmidt

function cm_eq_solution()
end

function cm_rand_params()
    spwp = rand(Uniform(0.20, 0.54)) #permanent wilting point
    sfc = spwp + 0.3                 #field capacity
    ep = 4.38                        #[mm/day] potential evaporation over land in mm/day, taken from [1]
    eo = 3.0                         #[mm/day] ocean evaporation rate
    ϵ = rand(Uniform(0.9, 1.1))      #numerical parameter from Rodriguez-Iturbe et al. (1991)
    r = rand(Uniform(1.9, 2.1))      #numerical parameter from Rodriguez-Iturbe et al. (1991)
    α = rand(Uniform(0.1,0.5))       #land fraction
    nZr = 100.0                      #[mm] reservoir depth/"field storage capacity of the soil" [mm] - NEEDS MORE RESEARCH
    a = rand(Uniform(11.4, 15.6))    #numerical parameter from Bretherton et al. (2004)
    b = rand(Uniform(0.522,0.603))   #numerical parameter from Bretherton et al. (2004)
    wsat = rand(Uniform(65.0, 80.0)) #[mm] saturation water vapour pass derived from plots in Bretherton et al.(2004) - NEEDS MORE RESEARCH
    u = rand(Uniform(5.0, 10.0)) * m2mm(1)/s2day(1) #[mm/day] wind speed
    L = rand(Uniform(500.0, 5000.0) * km2mm(1) #[mm] domain size
end
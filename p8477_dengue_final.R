#dengue final


# model function

SEIR_vector <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # host population
    dSH = (NH/TLH) - (SH * (cVH * (IV/NH) + (1/TLH)))
    dEH = SH * cVH * (IV/NH) - (EH * ((1/TIIT) + (1/TLH)))
    dIH = (EH/TIIT) - (IH * ((1/TID) + (1/TLH)))
    dRH = (IH/TID) - (RH/TLH)
    
    # vector population
    
    dSV = e - SV(cHV * ((IH + IH_visit)/NH) + (1/TLV))
    dEV = SV(cHV * ((IH + IH_visit)/NH)) - EV((1/TEIT) + (1/TLV))
    dIV = (EV/TEIT) - (IV/TLV)
    
    list(c(dSH, dEH, dIH, dRH, dSV, dEV, dIV))
    
  })
}

# parameter values for simulations


NH = 1e4    # host population
TLH = 600060  # host life span is 68.5 years or 600,060 hours
TIIT = 5  # intrinsic incubation period is 5 days
TEIT = 10   # extrinsic incubation period is 10 days
MPP = 2    # Number of mosquitoes per person, MPP = 2 for simulation 1-3
e = 5000    # emerging rate of adult mosquitoes, 5000 for simulation 1-3
IH_visit = 0  # visiting infectious host
TID = 3   # host infection duration is 3 days, AKA gamma
cVH = 0.75  # effective contact rate, vector to host is 0.75/day
cHV = 0.375   # effective contact rate, host to vector is 0.375/day

# run one of TLV for either wet or dry season

# parameter value for wet season 
TLV = 4   # vector life span is 4 days in wet season

# parameter value for dry season
TLV = 3   # vector life span is 3 days in dry season

# initial states (not directly from paper)

IH = 0  # initial number of infections in humans
EH = 1
RH = 
SH = NH - EH - IH


IV = 1  # initial number of infections in mosquitoes
SV = 

sim_season = ode()

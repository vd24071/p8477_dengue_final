#dengue final

library(deSolve)

# model function

SEIR_vector <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # host population
    dSH = (NH/TLH) - (SH * (cVH * (IV/NH) + (1/TLH)))
    dEH = (SH * cVH * (IV/NH)) - (EH * ((1/TIIT) + (1/TLH)))
    dIH = (EH/TIIT) - (IH * ((1/TID) + (1/TLH)))
    dRH = (IH/TID) - (RH/TLH)
    
    # vector population
    
    dSV = e - (SV * (cHV * ((IH + IH_visit)/NH) + (1/TLV)))
    dEV = (SV * (cHV * ((IH + IH_visit)/NH))) - (EV * ((1/TEIT) + (1/TLV)))
    dIV = (EV/TEIT) - (IV/TLV)
    
    list(c(dSH, dEH, dIH, dRH, dSV, dEV, dIV))
    
  })
}

# parameter values for simulations


NH = 1e4    # host population
TLH = 600060 # host life span is 68.5 years or 600,060 hours or 25002.5 days
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
TLV_wet = 4   # vector life span is 4 days in wet season

# parameter value for dry season
TLV_dry = 3   # vector life span is 3 days in dry season

# initial states (not directly from paper)

IH = 0  # initial number of infections in humans
EH = 1  # initial number of exposed in humans
RH = 0 # ???
SH = NH - EH - IH - RH

IV = 1  # initial number of infections in mosquitoes
EV = 1 # initial number of exposed in mosquitoes
SV = e # same as emerging rate???

times = 1:(365) # ???
  
state = c(SH = SH, EH = EH, IH = IH, RH = RH,
          SV = SV, EV = EV, IV = IV)

param_dry = c(TLH = TLH, TIIT = TIIT, TEIT = TEIT, MPP = MPP, e = e,
               IH_visit = IH_visit, TID = TID, cVH = cVH, cHV = cHV,
               TLV = TLV_dry)

sim_dry = ode(y=state,times=times,func=SEIR_vector,parms=param_dry)

matplot(sim_dry[,'time'], sim_dry[,c('SH','IH')], type = 'l', xlim = c(0,365), lwd = 1, col = c('blue', 'red'),
        lty = 1, main = 'Humans', cex.main = 1, ylab = 'Numbers', xlab = 'Time')
legend('topright',cex=1,seg.len = 1,
       legend=c('SH','IH'),
       lty=c(1,1),lwd=c(1,1),
       col=c('blue','red'),bty='n')

matplot(sim_dry[,'time'], sim_dry[,c('SV','IV')], type = 'l', xlim = c(0,365), lwd = 1, col = c('blue', 'red'),
        lty = 1, main = 'Mosquitoes', cex.main = 1, ylab = 'Numbers', xlab = 'Time')
legend('topright',cex=1,seg.len = 1,
       legend=c('SV','IV'),
       lty=c(1,1),lwd=c(1,1),
       col=c('blue','red'),bty='n')

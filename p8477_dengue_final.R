#dengue final

library(deSolve)

# basic model, no insecticide use

SEIR_basic <- function(t, state, parameters) {
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
    
    # cumulative incidence
    
    dcumInci = (1/TIIT) * EH
    
    list(c(dSH, dEH, dIH, dRH, dSV, dEV, dIV, dcumInci))
    
  })
}


# parameter values for simulations

NH = 1e4    # host population
TLH = 68.5*365 # host life span is 68.5 years or 600,060 hours or 25002.5 days
TIIT = 5  # intrinsic incubation period is 5 days
TEIT = 10   # extrinsic incubation period is 10 days
MPP = 2    # Number of mosquitoes per person, MPP = 2 for simulation 1-3
e = 5000    # emerging rate of adult mosquitoes, 5000 for simulation 1-3
IH_visit = 0  # visiting infectious host
TID = 3   # host infection duration is 3 days, AKA gamma
cVH = 0.75  # effective contact rate, vector to host is 0.75/day
cHV = 0.375   # effective contact rate, host to vector is 0.375/day

# average TLV

TLV = 4

# initial states (not directly from paper)

EH = 0  # initial number of exposed in humans 
IH = 0  # initial number of infections in humans 
RH = 0 # ???
SH = NH - EH - IH - RH

SV = 20000 # initial number of susceptible mosquitoes
EV = 0 # initial number of exposed in mosquitoes
IV = 1 # initial number of infectious in mosquitoes 

cumInci = 0 # initial cumulative incidence

times=seq(1,365) # ???

state = c(SH = SH, EH = EH, IH = IH, RH = RH,
          SV = SV, EV = EV, IV = IV, cumInci = cumInci)

param = c(TLH = TLH, TIIT = TIIT, TEIT = TEIT, MPP = MPP, e = e,
          IH_visit = IH_visit, TID = TID, cVH = cVH, cHV = cHV,
          TLV = TLV)

sim_basic = ode(y=state,times=times,func=SEIR_basic_no_fog,parms=param)

#matplot(sim_basic[,'time'], sim_basic[,'IV'], type = 'l', lwd = 1, col = 'blue', lty = 1)
matplot(sim_basic[,'time'], sim_basic[,'IH'], type = 'l', xlim = c(0,365), lwd = 1, col = 'red',
        lty = 1, main = 'Human Infections of Dengue', cex.main = 1, ylab = 'Prevalence of Dengue', xlab = 'Time')
legend('topright',cex=1,seg.len = 1,
       legend=c('IH'),
       lty=c(1,1),lwd=c(1,1),
       col=c('red'),bty='n')


# basic model with insecticide use

parms_no_fog = c(TLH = TLH, TIIT = TIIT, TEIT = TEIT, MPP = MPP, e = e,
                  IH_visit = IH_visit, TID = TID, cVH = cVH, cHV = cHV,
                  TLV = TLV)
state_no_fog = c(SH = SH, EH = EH, IH = IH, RH = RH,
          SV = SV, EV = EV, IV = IV, cumInci = cumInci)

times_no_fog=1;


parms_fog = c(TLH = TLH, TIIT = TIIT, TEIT = TEIT, MPP = MPP, e = e,
              IH_visit = IH_visit, TID = TID, cVH = cVH, cHV = cHV,
              TLV = TLV)

sim_no_fog = ode(y=state_no_fog,times=times_no_fog,func=SEIR_basic,parms=parms_no_fog);

res[5,]=tail(sim[,'cumInci'],1) # get daily cases

ts=seq(1,365,by=1)
res=matrix(0,length(ts),1)
prev=matrix()

for (i in 2:length(ts)){
  
  times_no_fog=seq(1,ts[i],by=1);
  times_fog=seq(ts[i],365,by=1);
  
  parms_fog = c(TLH = TLH, TIIT = TIIT, TEIT = TEIT, MPP = MPP, e = e,
                IH_visit = IH_visit, TID = TID, cVH = cVH, cHV = cHV,
                TLV = TLV)
  
  sim_no_fog = ode(y=state_no_fog,times=times_no_fog,func=SEIR_basic,parms=parms_no_fog);

  state_fog=c(SH=tail(sim_no_fog[,'SH'],1),EH=tail(sim_no_fog[,'EH'],1),IH=tail(sim_no_fog[,'IH'],1),
           RH=tail(sim_no_fog[,'RH'],1),SV=0.4*tail(sim_no_fog[,'SV'],1),EV=0.4*tail(sim_no_fog[,'EV'],1),
           IV=0.4*tail(sim_no_fog[,'IV'],1), cumInci=tail(sim_no_fog[,'cumInci'],1));
  
  sim_fog=ode(y=state_fog,times=times_fog,func=SEIR_basic,parms=parms_fog);
  
  sim=rbind(sim_no_fog,sim_fog[-1,])

  # save the result before exiting the loop
  res[i,]=tail(sim[,'cumInci'], 1) # get daily cases
}

#day 1 fogging
times_fog=seq(1,365,by=1);


#matplot(sim_basic[,'time'], sim_basic[,'IV'], type = 'l', lwd = 1, col = 'blue', lty = 1)
matplot(sim_basic_fog[,'time'], sim_basic_fog[,'IH'], type = 'l', xlim = c(0,365), lwd = 1, col = 'red',
        lty = 1, main = 'Human Infections of Dengue', cex.main = 1, ylab = 'Prevalence of Dengue', xlab = 'Time')
legend('topright',cex=1,seg.len = 1,
       legend=c('IH'),
       lty=c(1,1),lwd=c(1,1),
       col=c('red'),bty='n')

##################################################################

# model function for seasonality w/o insecticide use

SEIR_season_no_fog <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # real-time TLV based on season forcing (corrected) 
    #correction.factor = (1/365*((1+b.term)*150+(1-b.term)*215))
    #TLV = TLV/correction.factor * (1 + b.term * Term[time])
    
    # real-time TLV based on seasonality forcing (uncorrected) 
    TLV = TLV_avg + (TLV.term * Term[times])
    
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

## SETTING UP THE Season TERM-TIME FUNCTION:
# term-time forcing
wet_season=c(1:150);  # the days in a year of wet season
times=seq(1,365); # in day for 100 years
Term=rep(1,length(times));  # initialize a vector to store the Term
ind=(1:length(Term) %% 365) %in% wet_season  # find those days that are in wet season
Term[ind]=-1; # set them to -1

ind[1]

ind[160]

# parameter values for simulations


NH = 1e4    # host population
TLH = 68.5*365 # host life span is 68.5 years or 600,060 hours or 25002.5 days
TIIT = 5  # intrinsic incubation period is 5 days
TEIT = 10   # extrinsic incubation period is 10 days
MPP = 2    # Number of mosquitoes per person, MPP = 2 for simulation 1-3
e = 5000    # emerging rate of adult mosquitoes, 5000 for simulation 1-3
IH_visit = 0  # visiting infectious host
TID = 3   # host infection duration is 3 days, AKA gamma
cVH = 0.75  # effective contact rate, vector to host is 0.75/day
cHV = 0.375   # effective contact rate, host to vector is 0.375/day

TLV_avg = 3.5
TLV.term = 0.5

# run one of TLV for either wet or dry season

# parameter value for wet season 
#TLV_wet = 4   # vector life span is 4 days in wet season

# parameter value for dry season
#TLV_dry = 3   # vector life span is 3 days in dry season

# initial states (not directly from paper)

EH = 0  # initial number of exposed in humans (0.5*SH)
IH =  0 # initial number of infections in humans (0.25*EH)
RH = 0 # ???
SH = NH - EH - IH - RH

SV = 20000 # ratio of mosquitoes to humans
EV = 0 # initial number of exposed in mosquitoes (0.5*SV)
IV = 1 # initial number of infectious in mosquitoes 


state = c(SH = SH, EH = EH, IH = IH, RH = RH,
          SV = SV, EV = EV, IV = IV)

param = c(TLH = TLH, TIIT = TIIT, TEIT = TEIT, MPP = MPP, e = e,
               IH_visit = IH_visit, TID = TID, cVH = cVH, cHV = cHV,
               TLV = TLV)

sim_season = ode(y=state,times=times,func=SEIR_season_no_fog,parms=param)



# example attempts to plot model below

matplot(sim_season[,'time'], sim_season[,'IH'], type = 'l', xlim = c(0,365), lwd = 1, col = 'red',
        lty = 1, main = 'Human Infections of Dengue', cex.main = 1, ylab = 'Prevalence of Dengue', xlab = 'Time')
legend('topright',cex=1,seg.len = 1,
       legend='IH',
       lty=c(1,1),lwd=c(1,1),
       col='red',bty='n')

matplot(sim_dry[,'time'], sim_dry[,c('SH','EH', 'IH')], type = 'l', xlim = c(0,100), lwd = 1, col = c('green', 'blue', 'red'),
        lty = 1, main = 'Humans', cex.main = 1, ylab = 'Numbers', xlab = 'Time')
legend('topright',cex=1,seg.len = 1,
       legend=c('SH', 'EH', 'IH'),
       lty=c(1,1),lwd=c(1,1),
       col=c('green', 'blue','red'),bty='n')


matplot(sim_dry[,'time'], sim_dry[,c('SV','IV')], type = 'l', xlim = c(0,365), lwd = 1, col = c('blue', 'red'),
        lty = 1, main = 'Mosquitoes', cex.main = 1, ylab = 'Numbers', xlab = 'Time')
legend('topright',cex=1,seg.len = 1,
       legend=c('SV','IV'),
       lty=c(1,1),lwd=c(1,1),
       col=c('blue','red'),bty='n')


# wet season

param_wet = c(TLH = TLH, TIIT = TIIT, TEIT = TEIT, MPP = MPP, e = e,
              IH_visit = IH_visit, TID = TID, cVH = cVH, cHV = cHV,
              TLV = TLV_wet)

sim_wet = ode(y=state,times=times,func=SEIR_vector,parms=param_wet)

matplot(sim_wet[,'time'], sim_wet[,c('SH','EH', 'IH')], type = 'l', xlim = c(0,100), lwd = 1, col = c('green', 'blue', 'red'),
        lty = 1, main = 'Humans', cex.main = 1, ylab = 'Numbers', xlab = 'Time')
legend('topright',cex=1,seg.len = 1,
       legend=c('SH','EH','IH'),
       lty=c(1,1),lwd=c(1,1),
       col=c('green','blue','red'),bty='n')

matplot(sim_wet[,'time'], sim_wet[,c('SV','IV')], type = 'l', xlim = c(0,365), lwd = 1, col = c('blue', 'red'),
        lty = 1, main = 'Mosquitoes', cex.main = 1, ylab = 'Numbers', xlab = 'Time')
legend('topright',cex=1,seg.len = 1,
       legend=c('SV','IV'),
       lty=c(1,1),lwd=c(1,1),
       col=c('blue','red'),bty='n')


# need to simulate wet season for 5 months, then dry season for 7 months (rest of year)

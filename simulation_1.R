# Simulation 1

#load deSolve
library(deSolve)

#function for basic SEIR model

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

NH = 10000    # host population
TLH = 68.5*365 # host life span is 68.5 years or 600,060 hours or 25002.5 days
TIIT = 5  # intrinsic incubation period is 5 days
TEIT = 10   # extrinsic incubation period is 10 days
MPP = 2    # Number of mosquitoes per person, MPP = 2 for simulation 1-3
e = 5000    # emerging rate of adult mosquitoes, 5000 for simulation 1-3
IH_visit = 0  # visiting infectious host
TID = 3   # host infection duration is 3 days, AKA gamma
cVH = 0.75  # effective contact rate, vector to host is 0.75/day
cHV = 0.375   # effective contact rate, host to vector is 0.375/day

# TLV

TLV = 4

# initial states (from N&R paper)

EH = 0  # initial number of exposed in humans 
IH = 0  # initial number of infections in humans 
RH = 0 # initial number of recovered in humans
SH = NH - EH - IH - RH # initial number of suscepbile humans

SV = 20000 # initial number of susceptible mosquitoes
EV = 0 # initial number of exposed in mosquitoes
IV = 1 # initial number of infectious in mosquitoes 

cumInci = 0 # initial cumulative incidence

times=seq(1,365) # daily time step

state = c(SH = SH, EH = EH, IH = IH, RH = RH,
          SV = SV, EV = EV, IV = IV, cumInci = cumInci)

param = c(TLH = TLH, TIIT = TIIT, TEIT = TEIT, MPP = MPP, e = e,
          IH_visit = IH_visit, TID = TID, cVH = cVH, cHV = cHV,
          TLV = TLV)

sim_basic = ode(y=state,times=times,func=SEIR_basic,parms=param)

tail(sim_basic[,'cumInci'], 1)

sim_basic[which.max(sim_basic[,'IH'])]

# basic model with insecticide use

parms_no_fog = c(TLH = TLH, TIIT = TIIT, TEIT = TEIT, MPP = MPP, e = e,
                 IH_visit = IH_visit, TID = TID, cVH = cVH, cHV = cHV,
                 TLV = TLV)
state_no_fog = c(SH = SH, EH = EH, IH = IH, RH = RH,
                 SV = SV, EV = EV, IV = IV, cumInci = cumInci)

parms_fog = c(TLH = TLH, TIIT = TIIT, TEIT = TEIT, MPP = MPP, e = e,
              IH_visit = IH_visit, TID = TID, cVH = cVH, cHV = cHV,
              TLV = TLV)


ts=seq(1,365,by=1)
res=matrix(0,length(ts),2)
res[1:365,]=ts

for (i in 2:(length(ts)-1)){
  
  times_no_fog=seq(1,ts[i],by=1);
  times_fog=seq(ts[i],365,by=1);
  
  sim_no_fog = ode(y=state_no_fog,times=times_no_fog,func=SEIR_basic,parms=parms_no_fog);
  
  state_fog=c(SH=tail(sim_no_fog[,'SH'],1),EH=tail(sim_no_fog[,'EH'],1),IH=tail(sim_no_fog[,'IH'],1),
              RH=tail(sim_no_fog[,'RH'],1),SV=0.4*tail(sim_no_fog[,'SV'],1),EV=0.4*tail(sim_no_fog[,'EV'],1),
              IV=0.4*tail(sim_no_fog[,'IV'],1), cumInci=tail(sim_no_fog[,'cumInci'],1));
  
  sim_fog=ode(y=state_fog,times=times_fog,func=SEIR_basic,parms=parms_fog);
  
  sim=rbind(sim_no_fog,sim_fog[-1,])
  
  # save the result before exiting the loop
  res[i,2]=tail(sim[,'cumInci'], 1) # get cumulative incidence
}

#day 1 fogging
day_1_times_fog=seq(1,365,by=1);

state_fog=c(SH=SH,EH=EH,IH=IH,
            RH=RH,SV=0.4*SV,EV=0.4*EV,
            IV=0.4*IV, cumInci=cumInci);

day_1_sim_fog=ode(y=state_fog,times=day_1_times_fog,func=SEIR_basic,parms=parms_fog);

res[1,2]=tail(day_1_sim_fog[,'cumInci'],1) # get daily cases

#day 365 fogging

day_365_times_no_fog=seq(1,365,by=1);

day_365_sim_no_fog = ode(y=state_no_fog,times=day_365_times_no_fog,func=SEIR_basic,parms=parms_no_fog);

res[365,2]=tail(day_365_sim_no_fog[,'cumInci'], 1)

min(res[,2])
res[which.min(res[,2])]

View(res)


########################################

#Figure 1:

#Figure 1a:
# Day 169 fogging (optimal day)

times_no_fog=seq(1,169,by=1);
times_fog=seq(169,365,by=1);

sim_no_fog = ode(y=state_no_fog,times=times_no_fog,func=SEIR_basic,parms=parms_no_fog);

state_fog=c(SH=tail(sim_no_fog[,'SH'],1),EH=tail(sim_no_fog[,'EH'],1),IH=tail(sim_no_fog[,'IH'],1),
            RH=tail(sim_no_fog[,'RH'],1),SV=0.4*tail(sim_no_fog[,'SV'],1),EV=0.4*tail(sim_no_fog[,'EV'],1),
            IV=0.4*tail(sim_no_fog[,'IV'],1), cumInci=tail(sim_no_fog[,'cumInci'],1));

sim_fog=ode(y=state_fog,times=times_fog,func=SEIR_basic,parms=parms_fog);

sim_optimal=rbind(sim_no_fog,sim_fog[-1,])



plot(sim_basic[,'time'], sim_basic[,'IH'], type = 'l', lwd = 5, col = 'black',
        lty = 1, cex.axis = 1.3, cex.lab = 1.3, ylab = 'Prevalence of Dengue', xlab = 'Days', xaxs = 'i', xlim = c(0,365), ylim = c(0,300))
lines(sim_optimal[,'time'], sim_optimal[,'IH'], type = 'l', lwd = 5, lty = 3, col = 'red')
legend('right',cex=1.1,seg.len = 2,
       legend=c('No Fogging', 'Fogging on Day 169'),
       lty=c(1,3),lwd=c(5,5),
       col=c('black', 'red'),bty='n')
legend('topleft', cex = 1.4,
       legend = 'Wet Season',
       bty = 'n')

rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "slategray1", border = NA)
par(new = TRUE)


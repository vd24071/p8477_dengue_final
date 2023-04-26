# model function for seasonality 

#load deSolve
library(deSolve)

SEIR_season <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    #pull from seasonality vector
    
    TLV = seasonality[t]
    
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

# initial states 

EH = 0  # initial number of exposed in humans (0.5*SH)
IH =  0 # initial number of infections in humans (0.25*EH)
RH = 0 # initial number of recovered in humans
SH = NH - EH - IH - RH # initial number of suscepbile humans

SV = 20000 # initial number of susceptible mosquitoes
EV = 0 # initial number of exposed in mosquitoes
IV = 1 # initial number of infectious in mosquitoes 

cumInci = 0 # initial cumulative incidence

times=seq(1,365) 


## SETTING UP THE Season TERM-TIME FUNCTION for 4 months:
# term-time forcing

#run this for 4 months
seasonality = matrix(4,nrow=365,ncol=1)
seasonality[(30*4):365] = 3

#run this for 5 months
seasonality = matrix(4,nrow=365,ncol=1)
seasonality[(30*5):365] = 3

#run this for 6 months
seasonality = matrix(4,nrow=365,ncol=1)
seasonality[(30*6):365] = 3


state = c(SH = SH, EH = EH, IH = IH, RH = RH,
          SV = SV, EV = EV, IV = IV, cumInci = cumInci)

param = c(TLH = TLH, TIIT = TIIT, TEIT = TEIT, MPP = MPP, e = e,
          IH_visit = IH_visit, TID = TID, cVH = cVH, cHV = cHV)

sim_season = ode(y=state,times=times,func=SEIR_season,parms=param)


tail(sim_season[,'cumInci'], 1)

sim_season[which.max(sim_season[,'IH'])]


# seasonality with fogging


parms_no_fog = c(TLH = TLH, TIIT = TIIT, TEIT = TEIT, MPP = MPP, e = e,
                 IH_visit = IH_visit, TID = TID, cVH = cVH, cHV = cHV)

state_no_fog = c(SH = SH, EH = EH, IH = IH, RH = RH,
                 SV = SV, EV = EV, IV = IV, cumInci = cumInci)

parms_fog = c(TLH = TLH, TIIT = TIIT, TEIT = TEIT, MPP = MPP, e = e,
              IH_visit = IH_visit, TID = TID, cVH = cVH, cHV = cHV)


ts=seq(1,365,by=1)
res=matrix(0,length(ts),2)
res[1:365,]=ts

for (i in 2:(length(ts)-1)){
  
  times_no_fog=seq(1,ts[i],by=1);
  times_fog=seq(ts[i],365,by=1);
  
  sim_season_no_fog = ode(y=state_no_fog,times=times_no_fog,func=SEIR_season,parms=parms_no_fog);
  
  state_season_fog=c(SH=tail(sim_season_no_fog[,'SH'],1),EH=tail(sim_season_no_fog[,'EH'],1),IH=tail(sim_season_no_fog[,'IH'],1),
                     RH=tail(sim_season_no_fog[,'RH'],1),SV=0.4*tail(sim_season_no_fog[,'SV'],1),EV=0.4*tail(sim_season_no_fog[,'EV'],1),
                     IV=0.4*tail(sim_season_no_fog[,'IV'],1), cumInci=tail(sim_season_no_fog[,'cumInci'],1));
  
  sim_season_fog=ode(y=state_season_fog,times=times_fog,func=SEIR_season,parms=parms_fog);
  
  sim_seas=rbind(sim_season_no_fog,sim_season_fog[-1,])
  
  # save the result before exiting the loop
  res[i,2]=tail(sim_seas[,'cumInci'], 1) # get cumulative incidence
}

#day 1 fogging
day_1_times_fog=seq(1,365,by=1);

state_fog=c(SH=SH,EH=EH,IH=IH,
            RH=RH,SV=0.4*SV,EV=0.4*EV,
            IV=0.4*IV, cumInci=cumInci);

day_1_sim_seas_fog=ode(y=state_fog,times=day_1_times_fog,func=SEIR_season,parms=parms_fog);

res[1,2]=tail(day_1_sim_seas_fog[,'cumInci'],1) # get daily cases

#day 365 fogging

day_365_times_no_fog=seq(1,365,by=1);

day_365_sim_seas_no_fog = ode(y=state_no_fog,times=day_365_times_no_fog,func=SEIR_season,parms=parms_no_fog);

res[365,2]=tail(day_365_sim_seas_no_fog[,'cumInci'], 1)

min(res[,2])
res[which.min(res[,2])]

View(res)

#Figure 1b:

# Day 112 fogging (optimal day)

times_no_fog=seq(1,112, by = 1)
times_fog=seq(112,365,by=1);

sim_no_fog = ode(y=state_no_fog,times=times_no_fog,func=SEIR_season,parms=parms_no_fog);

state_fog=c(SH=tail(sim_no_fog[,'SH'],1),EH=tail(sim_no_fog[,'EH'],1),IH=tail(sim_no_fog[,'IH'],1),
            RH=tail(sim_no_fog[,'RH'],1),SV=0.4*tail(sim_no_fog[,'SV'],1),EV=0.4*tail(sim_no_fog[,'EV'],1),
            IV=0.4*tail(sim_no_fog[,'IV'],1), cumInci=tail(sim_no_fog[,'cumInci'],1));

sim_fog=ode(y=state_fog,times=times_fog,func=SEIR_season,parms=parms_fog);

sim_optimal_seas=rbind(sim_no_fog,sim_fog[-1,])

 

plot(sim_season[,'time'], sim_season[,'IH'], type = 'l',lwd = 5, col = 'black',
     lty = 1, cex.axis = 1.3, cex.lab = 1.3, ylab = 'Prevalence of Dengue', xlab = 'Days', xaxs = "i", xlim = c(0,365), ylim = c(0,300)) #  
lines(sim_optimal_seas[,'time'], sim_optimal_seas[,'IH'], type = 'l', lwd = 5, lty = 3, col = 'red')
legend('right',cex=1.1,seg.len = 2,
       legend=c('No Fogging', 'Fogging on Day 112'),
       lty=c(1,3),lwd=c(5,5),
       col=c('black', 'red'),bty='n')
legend('topleft', cex = 1.4,
       legend = 'Wet Season',
       bty = 'n')
legend('topright', cex = 1.4,
       legend = 'Dry Season',
       bty = 'n')

rect(0, par("usr")[3], 150, par("usr")[4], col = "slategray1", border = NA)
par(new = TRUE)

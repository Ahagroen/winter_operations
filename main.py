from dataclasses import dataclass
from turtle import clear
from typing import Literal
from gurobipy import GRB, LinExpr, Model, quicksum

model = Model("runway_winter")


# Knowns -> 

# Each runway becomes unsafe 1 time in the planning horizon
# Datasets ***
# Aircraft []
# Target time, 
# last possible time, 
# cost coefficient, 
# Minimum Sep Time between a and b (Same runway) - matrix
# Minimum sep time between a and cleaning
# Snow Removal Groups []
# Runways []
# Time to runway unsafety -120 min
# Time to remove snow - 20 min
# Setup time between runway clearing starts (including plow time and travel time) - 20min + travel time
#
#
# Decision Vars
# AC takeoff and landing times
# Start time of runway clearings
# If I happens before J (Binary)
# If ac a is on runway r (Binary)
# If ac a and ac b are on same runway (binary)
# If runway is cleared by snow removal group G (Binary)
# If runway a and b are cleared by the same group (Binary)


@dataclass
class Aircraft:
    target_time:int
    ac_class:Literal["Medium","Heavy","Super"]
    direction:Literal["Takeoff","Landing"]
    
    
#Inputs:
aircraft = []
planning_horizon = 0
runway_unsafe_times = [1500,1500]


#Constants
def last_time(ac:Aircraft):
    if ac.direction == "Takeoff":
        return planning_horizon
    else:
        return ac.target_time+20*60

def seperation_time(ac1:Aircraft,ac2:Aircraft):
    match ac1.direction:
        case "Landing":
            a = 0
        case "Takeoff":
            a = 3
    match ac1.ac_class:
        case "Medium":
            a += 0
        case "Heavy":
            a += 1
        case "Super":
            a += 2

    match ac2.direction:
        case "Landing":
            b = 0
        case "Takeoff":
            b = 3
    match ac2.ac_class:
        case "Medium":
            b += 0
        case "Heavy":
            b += 1
        case "Super":
            b += 2
    
    matrix = [[69,60,60,75,75,75],[157,96,96,75,75,75],[180,120,120,180,120,120],[60,60,60,60,60,60],[60,60,60,120,90,90],[180,120,120,180,120,120]]
    return matrix[a][b]


sep_a_cleaning = 60 #No idea, not in literature
t_snow_removal = 20*60
cost_coefficient = {"Medium":1,"Heavy":3,"Super":4}
runway_travel_matrix = [[0,600,1200],[900,0,1200],[1500,1500,0]] #This +20min is the value of Q

#Run-specific Constants
snow_removal_groups:int = 0
runways:int = 0


event_times = {}

for i in aircraft:
    event_times[aircraft.index(i)] = model.addVar(name="x_"+str(aircraft.index(i)))

clearing_times = {}
for r in range(0,runways):
    clearing_times[r] = model.addVar(name="v"+str(r))

count = 0
i_before_j = {}
for i in aircraft + list(range(runways)):
    i_before_j[count] = model.addVar(vtype=GRB.BINARY,name="delta_"+str(count))
    count += 1

yar = {}
for i in aircraft:
    for j in range(runways):
        yar[aircraft.index(i),j] = model.addVar(vtype=GRB.BINARY,name="y_"+str(aircraft.index(i))+"_"+str(j))
zij = {}
for i in aircraft:
    for j in aircraft:
        if j != i:
            zij[aircraft.index(i),aircraft.index(j)] = model.addVar(vtype=GRB.BINARY,name="z_"+str(aircraft.index(i))+str(aircraft.index(j)))

rho = {}
for i in range(runways):
    for j in range(snow_removal_groups):
        rho[i,j] = model.addVar(vtype=GRB.BINARY,name="rho_"+str(i)+"_"+str(j))

psi = {}
for i in range(runways):
    for j in range(runways):
        if i != j:
            psi[i,j] = model.addVar(vtype=GRB.BINARY,name="psi_"+str(i)+"_"+str(j))


obj = LinExpr()
for i in aircraft:
    obj+= cost_coefficient[i.ac_class]*(event_times[aircraft.index(i)] - i.target_time)

model.setObjective(obj,GRB.MINIMIZE)



for i in aircraft:
    model.addConstr(i.target_time < event_times[aircraft.index(i)] < last_time(i),"C0_"+str(aircraft.index(i)))

for i in i_before_j.keys():
    for j in i_before_j.keys():
        if i != j:
            model.addConstr(i_before_j[i]+i_before_j[j] ==1,"C1_"+str(i)+"_"+str(j))

for i in aircraft:
    model.addConstr(quicksum(yar[aircraft.index(i),j] for j in range(runways)) == 1)

for i in aircraft:
    for j in aircraft:
        if i != j:
            for k in range(runways):
                model.addConstr(zij[aircraft.index(i),aircraft.index(j)] >= yar[aircraft.index(i),k] + yar[aircraft.index(j),k] -1)


for i in aircraft:
    for j in aircraft:
        if i != j:
            model.addConstr(event_times[aircraft.index(j)] >= event_times[aircraft.index(i)] + seperation_time(i,j)*zij[aircraft.index(i),aircraft.index(j)] - 100000*i_before_j[aircraft.index(j),aircraft.index(i)])

for i in range(runways):
    model.addConstr(quicksum(rho[i,j] for j in range(snow_removal_groups)) == 1)

for i in range(runways):
    for j in range(runways):
        if i != j:
            for k in range(snow_removal_groups):
                model.addConstr(psi[i,j] >= rho[i,k] + rho[j,k] -1 )


for i in range(runways):
    for j in range(runways):
        if i != j:
            model.addConstr(clearing_times[j] >= clearing_times[i] + (t_snow_removal+runway_travel_matrix[i][j])*psi[i,j]-100000*i_before_j[j,i])


for i in range(runways):
    for j in aircraft:
        model.addConstr(clearing_times[i] >= event_times[aircraft.index(j)] + sep_a_cleaning*yar[i,aircraft.index(j)]-100000*i_before_j[i,aircraft.index(j)])

for i in range(runways):
    for j in aircraft:
        model.addConstr(event_times[aircraft.index(j)] >= clearing_times[i] + t_snow_removal*yar[i,aircraft.index(j)]-100000*i_before_j[aircraft.index(j),i])


for i in range(len(aircraft)):
    for j in range(runways):
        model.addConstr(event_times[i] <= runway_unsafe_times[j] + 100000*(1-yar[i,j]) + 100000*i_before_j[j,i])




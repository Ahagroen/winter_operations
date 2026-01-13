from dataclasses import dataclass, field
from itertools import count
from typing import Literal
from gurobipy import GRB, LinExpr, Model, quicksum
import csv

from matplotlib import pyplot

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
    identifier: int = field(default_factory=count().__next__, init=False)
    target_time:int
    ac_class:Literal["Medium","Heavy","Super"]
    direction:Literal["Takeoff","Landing"]
    
def map_type(type):
    match type:
        case "H":
            return "Heavy"
        case "M":
            return "Medium"
        case "S":
            return "Super"
        case _:
            raise RuntimeError("Invalid type in csv")

def map_direction(direction):
    match direction:
        case "L":
            return "Landing"
        case "T":
            return "Takeoff"
        case _:
            raise RuntimeError("Invalid direction in csv")
def load_aircraft(filename:str)->list[Aircraft]:
    carry = []
    with open(filename,"r") as fs:
        reader = csv.reader(fs)
        for row in reader:
            carry.append(Aircraft(int(row[0])*60,map_type(row[1]),map_direction(row[2])))
    return carry
#Inputs:
planning_horizon = 140*60
runway_unsafe_times = [3000,1500]
aircraft = load_aircraft("schipol1d.csv")

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


sep_a_cleaning = 120 #No idea, not in literature
t_snow_removal = 20*60
cost_coefficient = {"Medium":1,"Heavy":3,"Super":4}
runway_travel_matrix = [[0,600,1200],[900,0,1200],[1500,1500,0]] #This +20min is the value of Q

#Run-specific Constants
snow_removal_groups:int = 1
runways:list[str] = ["A","B"]


event_times = {}

for i in aircraft:
    event_times[aircraft.index(i)] = model.addVar(name="x_"+str(aircraft.index(i)))

clearing_times = {}
for r in runways:
    clearing_times[r] = model.addVar(name="v"+str(r))

i_before_j = {}
for i in list(range(len(aircraft))) + runways:
    for j in list(range(len(aircraft))) + runways:
        if i != j:
            i_before_j[i,j] = model.addVar(vtype=GRB.BINARY,name="delta_"+str(i)+"_"+str(j)) 

yar = {}
for i in aircraft:
    for j in runways:
        yar[aircraft.index(i),j] = model.addVar(vtype=GRB.BINARY,name="y_"+str(aircraft.index(i))+"_"+str(j))
zij = {}
for i in aircraft:
    for j in aircraft:
        if j != i:
            zij[aircraft.index(i),aircraft.index(j)] = model.addVar(vtype=GRB.BINARY,name="z_"+str(aircraft.index(i))+str(aircraft.index(j)))

rho = {}
for i in runways:
    for j in range(snow_removal_groups):
        rho[i,j] = model.addVar(vtype=GRB.BINARY,name="rho_"+str(i)+"_"+str(j))

psi = {}
for i in runways:
    for j in runways:
        if i != j:
            psi[i,j] = model.addVar(vtype=GRB.BINARY,name="psi_"+str(i)+"_"+str(j))


obj = LinExpr()
for i in aircraft:
    obj+= cost_coefficient[i.ac_class]*(event_times[aircraft.index(i)] - i.target_time)

model.setObjective(obj,GRB.MINIMIZE)



for i in aircraft:
    model.addConstr(i.target_time <= event_times[aircraft.index(i)])

for i in aircraft:
    model.addConstr(event_times[aircraft.index(i)] <= last_time(i),"C0_"+str(aircraft.index(i)))

for i in list(range(len(aircraft))) + runways:
    for j in list(range(len(aircraft))) + runways:
        if i != j:
            model.addConstr(i_before_j[i,j]+i_before_j[j,i] ==1,"C1_"+str(i)+"_"+str(j))

for i in aircraft:
    model.addConstr(quicksum(yar[aircraft.index(i),j] for j in runways) == 1)

for i in aircraft:
    for j in aircraft:
        if i != j:
            for k in runways:
                model.addConstr(zij[aircraft.index(i),aircraft.index(j)] >= yar[aircraft.index(i),k] + yar[aircraft.index(j),k] -1)


for i in aircraft:
    for j in aircraft:
        if i != j:
            model.addConstr(event_times[aircraft.index(j)] >= event_times[aircraft.index(i)] + seperation_time(i,j)*zij[aircraft.index(i),aircraft.index(j)] - 100000*i_before_j[aircraft.index(j),aircraft.index(i)])

for i in runways:
    model.addConstr(quicksum(rho[i,j] for j in range(snow_removal_groups)) == 1)

for i in runways:
    for j in runways:
        if i != j:
            for k in range(snow_removal_groups):
                model.addConstr(psi[i,j] >= rho[i,k] + rho[j,k] -1 )


for i in runways:
    for j in runways:
        if i != j:
            model.addConstr(clearing_times[j] >= clearing_times[i] + (t_snow_removal+runway_travel_matrix[runways.index(i)][runways.index(j)])*psi[i,j]-100000*i_before_j[j,i])


for i in runways:
    for j in aircraft:
        model.addConstr(clearing_times[i] >= event_times[aircraft.index(j)] + sep_a_cleaning*yar[aircraft.index(j),i]-100000*i_before_j[i,aircraft.index(j)])

for i in runways:
    for j in aircraft:
        model.addConstr(event_times[aircraft.index(j)] >= clearing_times[i] + t_snow_removal*yar[aircraft.index(j),i]-100000*i_before_j[aircraft.index(j),i])


for i in range(len(aircraft)):
    for j in runways:
        model.addConstr(event_times[i] <= runway_unsafe_times[runways.index(j)] + 100000*(1-yar[i,j]) + 100000*i_before_j[j,i])


model.optimize()

def match_runway(flight):
    return [x for x in runways if yar[flight.identifier,x].X == 1][0]

runway_flights = {}
delays = {}
for i in runways:
    runway_flights[i] = {"Takeoff":[],"Landing":[]}
for i in event_times.keys():
    runway_flights[match_runway(aircraft[i])][aircraft[i].direction].append(event_times[i].X)
    if event_times[i].X > aircraft[i].target_time:
        delays[aircraft[i].identifier] = (aircraft[i].target_time,event_times[i].X)
        print(f"Flight: {aircraft[i].identifier} occurs at {event_times[i].X} on runway {match_runway(aircraft[i])} - Delay of {event_times[i].X - aircraft[i].target_time}")
    else:
        print(f"Flight: {aircraft[i].identifier} occurs at {event_times[i].X} on runway {match_runway(aircraft[i])}")


counter = 0
for i in runway_flights.keys():
    pyplot.scatter(runway_flights[i]["Takeoff"],[counter]*len(runway_flights[i]["Takeoff"]),label="Takeoff")
    pyplot.scatter(runway_flights[i]["Landing"],[counter]*len(runway_flights[i]["Landing"]),label="Landing")
    pyplot.plot([clearing_times[i].X,clearing_times[i].X+20*60],[counter]*2,label="Clearing")
    counter +=1

pyplot.legend()
pyplot.show()

counter = 0
for i in delays.keys():
    pyplot.plot(range(int(delays[i][0]),int(delays[i][1])),range(0,int(delays[i][1]-delays[i][0])),label=i)

pyplot.legend()
pyplot.show()
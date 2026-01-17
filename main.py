from collections import Counter
from dataclasses import dataclass, field
from enum import Enum
from functools import cache
from itertools import count
from logging import DEBUG, INFO
from sys import stdout
from typing import Literal
from gurobipy import GRB, LinExpr, Model, quicksum
import csv
from loguru import logger
import time
from os import path

from matplotlib import pyplot



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
        
@cache
def load_aircraft(filename:str)->list[Aircraft]:
    carry = []
    with open(filename,"r") as fs:
        reader = csv.reader(fs)
        for row in reader:
            carry.append(Aircraft(int(row[0]),map_type(row[1]),map_direction(row[2])))
    return carry




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

class Mode(Enum):
    large_scale = 0
    delayed_arrival = 1
    infeasible_landing_delay = 2
    delayed_after_plow = 4
    too_short_one_plow = 5
    validation = 6


def callback(model, where):
    if where == GRB.Callback.MIPNODE:
        # Get model objective
        obj = model.cbGet(GRB.Callback.MIPNODE_OBJBST)

        # Has objective changed?
        if abs(obj - model._cur_obj) > 1e-8:
            # If so, update incumbent and time
            model._cur_obj = obj
            model._time = time.time()

    # Terminate if objective has not improved in 20s
    if time.time() - model._time > 30:
        model.terminate()



def run_model(mode:Mode,returns:bool,plowing_time:int=20*60,num_runways:int=3,num_plows:int=2,event_margin:int=20*60,sep_fuzz_factor:float=0,runway_unsafe_times:list[int]=[1500,1500,1500,1500,1500],file_name:str=""):
    model = Model("runway_winter")
    model.setParam('TimeLimit', 5*60)
    # model.setParam("OutputFlag",0)
    
    #set up early stopping
    match mode:
        case Mode.large_scale:
            #Inputs - Large Scale:
            planning_horizon = 120*60
            aircraft = load_aircraft("schipol1d.csv")
        case Mode.delayed_arrival:
            #Inputs - Delayed arrival:
            planning_horizon = 60*60
            runway_unsafe_times = [1200,1500]
            num_plows = 1
            num_runways = 2
            aircraft = load_aircraft("test_file_delay.csv")
        case Mode.infeasible_landing_delay:
            #Inputs - Infeasible (too long landing delay):
            planning_horizon = 60*60
            runway_unsafe_times = [120,120]
            num_plows = 1
            num_runways = 2
            aircraft = load_aircraft("test_file_infeasible.csv")
        case Mode.delayed_after_plow:
            #Inputs - Delayed until after plowing:
            planning_horizon = 60*60
            runway_unsafe_times = [120,120]
            num_plows = 1
            num_runways = 2
            aircraft = load_aircraft("test_file_plow.csv")
        case Mode.too_short_one_plow:
            planning_horizon = 25*60
            runway_unsafe_times = [120,120]
            num_plows = 1
            num_runways = 2
            aircraft = load_aircraft("test_file_plow.csv")
        case Mode.validation:
            planning_horizon = 3600
            num_plows = 1
            num_runways = 2
            runway_unsafe_times = [600,1500]
            plowing_time = 1200
            aircraft = load_aircraft("munichc9.csv")
        case _:
            raise RuntimeError("Not Implemented")

    #Constants
    sep_a_cleaning = 120 #No idea, not in literature
    t_snow_removal = plowing_time
    cost_coefficient = {"Medium":1,"Heavy":3,"Super":4}
    runway_travel_matrix = [[0,600,1200,900,1200],[900,0,1200,1200,600],[1500,1500,0,600,900],[1500,900,900,0,1500],[600,1200,900,1200,0]] #This +20min is the value of Q

    #Run-specific Constants
    snow_removal_groups:int = num_plows
    runways:list[str] = ["A","B","C","D","E"]
    runways = runways[0:num_runways]

    def last_time(ac:Aircraft):
        return ac.target_time+event_margin

    #Decision Vars  
    event_times = {} #Xa
    for i in aircraft:
        event_times[aircraft.index(i)] = model.addVar(vtype=GRB.INTEGER,name="x_"+str(aircraft.index(i)))

    clearing_times = {} #Vr
    for r in runways:
        clearing_times[r] = model.addVar(vtype=GRB.INTEGER,name="v"+str(r))

    i_before_j = {}
    for i in list(range(len(aircraft))) + runways:
        for j in list(range(len(aircraft))) + runways:
            if i != j:
                i_before_j[i,j] = model.addVar(vtype=GRB.BINARY,name="delta_"+str(i)+"_"+str(j)) 
                #Pre-solve, based on paper
                if isinstance(i,int) and isinstance(j,int):
                    if aircraft[j].target_time > 20*60+aircraft[i].target_time:
                        i_before_j[i,j].start = 1
                    elif aircraft[i].target_time < aircraft[j].target_time and aircraft[i].ac_class == aircraft[j].ac_class:
                        i_before_j[i,j].start = 1

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
                model.addConstr(event_times[aircraft.index(j)] >= event_times[aircraft.index(i)] + (1+sep_fuzz_factor)*seperation_time(i,j)*zij[aircraft.index(i),aircraft.index(j)] - 100000*i_before_j[aircraft.index(j),aircraft.index(i)])

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

    model._cur_obj = float('inf')
    model._time = time.time()
    model.optimize(callback=callback)

    def match_runway(flight):
        return [x for x in runways if yar[flight.identifier,x].X == 1][0]

    runway_flights = {}
    delays = {}
    total_delay = 0
    for i in runways:
        runway_flights = {"Takeoff":[],"Landing":[]}
    for i in event_times.keys():
        runway_flights[aircraft[i].direction].append((event_times[i].X,runways.index(match_runway(aircraft[i]))))
        if event_times[i].X > aircraft[i].target_time:
            total_delay += event_times[i].X - aircraft[i].target_time
            delays[aircraft[i].identifier] = (aircraft[i].target_time,event_times[i].X)
            logger.debug(f"Flight: {aircraft[i].identifier} occurs at {event_times[i].X} on runway {match_runway(aircraft[i])} - Delay of {event_times[i].X - aircraft[i].target_time}")
        else:
            logger.debug(f"Flight: {aircraft[i].identifier} occurs at {event_times[i].X} on runway {match_runway(aircraft[i])}")
    logger.debug(f"Total (non cost adjusted) delay: {total_delay}, Cost adjusted delay: {model.ObjVal}")
    for i in runway_flights.keys():
        if len(runway_flights[i]) == 0:
            continue
        x,y = zip(*runway_flights[i])
        pyplot.scatter(x,y,label=i)
    for index,val in enumerate(runways):
        if clearing_times[val].X < planning_horizon:
            pyplot.plot([clearing_times[val].X,clearing_times[val].X+plowing_time],[index]*2,label="Clearing" if index == 0 else "",color="orange")
    pyplot.xlabel("Time",fontsize=14)
    pyplot.ylabel("Runway",fontsize=14)
    pyplot.title("Event Times on each runway",fontsize=20)
    pyplot.yticks(range(len(runways)),runways)
    pyplot.legend()
    if not returns:
        pyplot.show()
    else:
        pyplot.gcf().set_size_inches(12,4)
        pyplot.savefig(path.join("outputs",f"{num_runways}-{num_plows}-{file_name}-times.svg"),dpi=100,bbox_inches='tight')
        pyplot.close()
    carry = []
    for i in delays.keys():
        carry+=range(int(delays[i][0]),int(delays[i][1]))
        if delays[i][1]-delays[i][0] > 0:
            pyplot.bar(i,delays[i][1]-delays[i][0],label=i)
    pyplot.xticks(range(0,len(aircraft)))
    pyplot.xlabel("Aircraft Number",fontsize=14)
    pyplot.ylabel("Delay (s)",fontsize=14)
    pyplot.title("Flight Delay in Seconds",fontsize=20)
    if not returns:
        pyplot.show()
    else:
        pyplot.gcf().set_size_inches(16,6)
        pyplot.savefig(path.join("outputs",f"{num_runways}-{num_plows}-{file_name}-delays.svg"),dpi=100,bbox_inches='tight')
        pyplot.close()
    c = Counter(carry)
    for i in c.keys():
        pyplot.bar(i,c[i],color="blue")
    pyplot.xlabel("Time (s)",fontsize=14)
    pyplot.ylabel("Aircraft Waiting",fontsize=14)
    pyplot.title("Total Aircraft Waiting",fontsize=20)
    if not returns:
        pyplot.show()
    else:
        pyplot.gcf().set_size_inches(12,8)
        pyplot.savefig(path.join("outputs",f"{num_runways}-{num_plows}-{file_name}-waiting.svg"),dpi=100,bbox_inches='tight')
        pyplot.close()
        return total_delay,model.ObjVal

def discrete_sensitivity_analysis():
    results_runways = []
    snow_cond = [[1500,1500,1500,1500,1500],[1500,2400,3600,2400,1500],[0,1500,600,1500,0]]
    for i in [2,3,4,5]: #[2,3,4,5]
        for j in range(1,i):
            carry = []
            for k in range(0,3):
                logger.info(f"snow type run - {i} runways, {j} plows, snow condition {k}")
                carry.append((i,j,k,run_model(Mode.large_scale,True,num_runways=i,num_plows=j,runway_unsafe_times=snow_cond[k],file_name="snow-"+str(k))))
            results_runways.append(carry)
    for i in results_runways:
        logger.info(f"With {i[0][0]} runways and {i[0][1]} plows in snow condition {i[0][2]}: total non-weighted delay: {i[0][3][0]}, total weighted delay: {i[0][3][1]}")
        logger.info(f"With {i[0][0]} runways and {i[0][1]} plows in snow condition {i[1][2]}: total non-weighted delay: {i[1][3][0]}, total weighted delay: {i[1][3][1]}")
        logger.info(f"With {i[0][0]} runways and {i[0][1]} plows in snow condition {i[2][2]}: total non-weighted delay: {i[2][3][0]}, total weighted delay: {i[2][3][1]}")
        pyplot.plot([0,1,2],[i[0][3][1],i[1][3][1],i[2][3][1]],label=f"{i[0][0]} Runways - {i[0][1]} Groups")
    pyplot.xlabel("Snow Condition")
    pyplot.ylabel("Weighted Sum of Delay")
    pyplot.xticks(range(0,3),["Starting Snow","Snowing","Ending Snow"])
    pyplot.title("Sensitivity Analysis on Snow Scenarios")
    pyplot.legend()
    pyplot.show()
    

def plow_time_sensitivity_analysis():
    carry = []
    for i in range(10,32,2):
        logger.info(f"plow time run - {i} minutes")
        carry.append((i,run_model(Mode.large_scale,True,num_runways=4,num_plows=3,plowing_time=i*60,file_name="plow-time"+str(i))))
    for i in carry:
        logger.info(f"With plow-time {i[0]}, total non-weighted delay: {i[1][0]}, total weighted delay: {i[1][1]}")
    pyplot.plot(range(10,32,2),[x[1][1] for x in carry])
    pyplot.xlabel("Time to Plow a runway (minutes)")
    pyplot.ylabel("Total Weighted Delay (s)")
    pyplot.title("Sensitivity of Runway Plow time on Delay")
    pyplot.show()


def landing_space_sensitivity_analysis():
    carry = []
    for i in range(0,12,2):
        logger.info(f"sep time run - {i} percent extra")
        carry.append((i/10,run_model(Mode.large_scale,True,num_runways=4,num_plows=3,sep_fuzz_factor=i/10,file_name="sep-time"+str(100+i))))
    for i in carry:
        logger.info(f"With sep-time margin {i[0]}, total non-weighted delay: {i[1][0]}, total weighted delay: {i[1][1]}")
    pyplot.plot([x for x in list(range(100,220,20))],[x[1][1] for x in carry])
    pyplot.xlabel("Percentage of Clear Weather seperation time")
    pyplot.ylabel("Total Weighted Delay (s)")
    pyplot.title("Sensitivity of Sab on Delay")
    pyplot.show()


def landing_time_sensitivity_analysis():
    carry = []
    for i in range(10,30,2):
        logger.info(f"event margin run - {i} minutes")
        carry.append((i,run_model(Mode.large_scale,True,num_runways=4,num_plows=3,event_margin=i*60,file_name="event_margin"+str(i))))
    for i in carry:
        logger.info(f"With event margin {i[0]}, total non-weighted delay: {i[1][0]}, total weighted delay: {i[1][1]}")
    pyplot.plot([x for x in list(range(100,220,20))],[x[1][1] for x in carry])
    pyplot.xlabel("Allowable delay past scheduled (m)")
    pyplot.ylabel("Total Weighted Delay (s)")
    pyplot.title("Sensitivity of Allowable Delay on Delay")
    pyplot.show()


if __name__ == "__main__":
    logger.remove()
    logger.add(stdout,level=INFO)
    logger.add("runlog.log",level=DEBUG)

    #Sensitivity Analysis
    # landing_space_sensitivity_analysis()
    # plow_time_sensitivity_analysis()
    # discrete_sensitivity_analysis()
    # landing_space_sensitivity_analysis()

    #Verification
    # run_model(Mode.delayed_arrival,True,num_plows=1,num_runways=2,file_name="verification-delayed-arrival")
    # run_model(Mode.delayed_after_plow,True,num_plows=1,num_runways=2,file_name="verification-delayed-afterplow")
    # run_model(Mode.infeasible_landing_delay,True,num_plows=1,num_runways=2,file_name="verification-infeasible-delay")
    # run_model(Mode.too_short_one_plow,True,num_plows=1,num_runways=2,file_name="verification-one-plow")
    # run_model(Mode.large_scale,True,num_plows=1,runway_unsafe_times=[0,1500,600,1500,0],file_name="snow-0")
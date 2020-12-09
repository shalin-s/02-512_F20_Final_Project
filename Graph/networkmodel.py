import random
#im actually using my javascript code for the plots/everything, but here is basically what I did in python in case its of interest
#links represented as array of sets
#nodes given as an dict with (id, state)
S = 100
I = 0
R =0
D =0
infectionRate = 0.6
recoverRate = 0.4
nodes = {}
links = []
def simStep(nodes, links):
    #we iterate through the nodes and pull all nodes with "I" state
    infect = []
    for a in nodes:
        if nodes[a]=="I":
            infect.append(a)
    #we find all links with attached infected nodes
    linked = []
    for s in links:
        for nodeId in infect:
            if nodeId in s:
                linked+=list(s)

    #we infect neighboring nodes
    for i in linked:
        #we randomly choose susceptible/infected nodes to change states from our linked
        if nodes[i] == "S":
            if(random.random() >infectionRate):
                nodes[i] = "I"
                S-=1
                I+=1
    for i in infect:
        if(random.random() > recoverRate):
            if(random.random() < 0.05):
                nodes[i]="D"
                I-=1
                D+=1
            else:
                nodes[i]= "R"
                I-=1
                R+=1

while(I>0):
    simStep(nodes, links)


#Importing relevant librairies
import numpy as np
import random
import pathlib
import matplotlib.pyplot as plt
import matplotlib.patches as ptc
import matplotlib.animation as anim
import matplotlib.collections as collec

path = str(pathlib.Path(__file__).parent.resolve())

#For the system we consider a subgraph of Z2 of the following form :     
#                
#       _ _ _ _ 
#      |_|_|_|_|
#      |_|_|_|_|
#      |_|_|_|_|
#      |_|_|_|_|
#      <------->
#         m
#
#If a graph is named simply G it means the values are 0 or 1 depending on if the edges belongs to the graph, if it's named Ga it means the value are the radius of the tubes, if it's named Gn it encodes the quanties of softening agent (the absolute quantity not the concentration)

def CreateEmptyGraph(m,n): #takes the parameter m and n the width and the height of the graph and creates a graph of such dimension with value 0 everywhere
    G = []
    for i in range(2*n+1): #because of the shape of such graphs the stucture is a bit unsual 
        if i%2 == 0:
            G.append(np.array([0. for _ in range(m)]))  #If i is even this means the row is enconding horizontal edges they are m of them 
        else:
            G.append(np.array([0. for _ in range(m+1)])) #If i is odd we encode vertical edges meaning there is 1 more : m+1
    return G                                             #In total we need 2n+1 list of edges to encode everything           


def GenerateRandomGraph(m,n,p):   # To get a graph we proceed in two steps first we create a random of graph of size m,n meaning we only keep an edge with probability p
    G = []
    for i in range(2*n+1):
        if i%2 == 0:
            G.append(np.array([0 if i == 0 else int(random.random()<p) for j in range(m)]))  
        else:
            G.append(np.array([1 if (j==m//2 and i == 1) else (0 if j == 0 or j == m else int(random.random()<p)) for j in range(m+1)]))
    return G     #Some edges are always discared (the ones at the bottom and at the left and right limits) and the central bottom edges is always kept such that the graphs always have similar sturcures
   

def ClearUnconnectedBits(G):   # The second step is to remove any unconnected bits of the graphs from the main part
    n,m = len(G)//2,len(G[0])
    B = CreateEmptyGraph(m,n)   #In order to do that we do a depth for search (DFS), meaning we explore the graph from neighbour to neighbour from the origin marking edges as we goes
    B[1][m//2] = 1
    Stack = [(1,m//2)]
    while len(Stack) > 0:
        i,j = Stack.pop(0)
        P1,P2 = PointFromGraphIndice((i,j))    #Those two function are defined and explained below
        Neigh = NeighbourVertice(G,P1,P2)      #
        for (u,v) in Neigh:
            if B[u][v] == 0:
                B[u][v] = 1
                Stack.append((u,v))
    for i in range(len(G)):                 #at the end we remove all the unmarked edg
        for j in range(len(G[i])):
            G[i][j] = G[i][j]*B[i][j]
    return G

#A few useful functions to work on those type of graphs :

def VerticeIndice(G,P1,P2):    #take two point P1 and P2 and return the indices i,j that give this edge in the graph G, if the points are not connected or are outisde of the graph it return "Not Valid" 
    d = P1 - P2 
    n,m = len(G)//2,len(G[0])
    if (np.array_equal(d,np.array([0,1])) or np.array_equal(d,np.array([0,-1])) or np.array_equal(d,np.array([1,0])) or np.array_equal(d,np.array([-1,0]))) and 0<=P1[0]<=m and 0<=P2[0]<=m and 0<=P1[1]<=n and 0<=P2[1]<=n:
        if d[0] == 0:
            return (2*min(P1[1],P2[1])+1,P1[0])
        else:
            return (2*P1[1],min(P1[0],P2[0]))
    else:
        return "Not Valid"
    
def PointFromGraphIndice(I):    #Inverse of the previous function, we give indices in the graph and get the two points forming that edge       
    u,v = I
    if u%2 == 0:
        return np.array([v,u//2]),np.array([v+1,u//2])
    else:
        return np.array([v,u//2]),np.array([v,u//2+1])
        

def neighbour(G,i,j):       #return the neighbours of the vertice situated at i,j in the graph
    L = [np.array([i,j])]
    n,m = len(G)//2,len(G[0])
    if i + 1 < m:                                                   #for each of the 4 potentials neighbour we check wether it's oustide the graph or not, and then check if the edges connecting the points exists in G and add it to the list depending on that.
        u,v = VerticeIndice(G,np.array([i,j]),np.array([i+1,j]))    # the point itself is always counted as a neighbour
        if G[u][v]:
            L.append(np.array([i+1,j]))
    if 0 <= i-1:
        u,v = VerticeIndice(G,np.array([i-1,j]),np.array([i,j]))
        if G[u][v]:
            L.append(np.array([i-1,j]))
    if j + 1 < n:
        u,v = VerticeIndice(G,np.array([i,j]),np.array([i,j+1]))
        if G[u][v]:
            L.append(np.array([i,j+1]))
    if 0 <= j-1:
        u,v = VerticeIndice(G,np.array([i,j-1]),np.array([i,j]))
        if G[u][v]:
            L.append(np.array([i,j-1]))
    return L

def NeighbourVertice(G,P1,P2):       #Here we gives the neighbour of an edge given by two point P1 and P2, an edge's neighbour is all the edges with which it shares a point
    L = [VerticeIndice(G,P1,P2)]
    u,v = P1
    for P in neighbour(G,u,v):                   #We loop throught the neighbours of each points and add the adge if it's correct (note that the edges are added in there incides form)
        if not(np.array_equal(P1,P)):            #Same as for the vertices the edge itself is considered one of it's neighbour
            L.append(VerticeIndice(G,P1,P))
    a,b = P2
    for P in neighbour(G,a,b):
        if not(np.array_equal(P2,P)):
            L.append(VerticeIndice(G,P2,P))
    return L

def GenerateIniGa(G,a0):   #We create an inital graph Ga of the same shape as G putting the value a0 everywhere 
    Ga = []
    for i in range(len(G)):
        Ga.append(np.array([a0 if G[i][j] else 0 for j in range(len(G[i]))]))
    return Ga   

alpha = 10  #This parameter controls the stenght of the effect the larger it is the more the tube will grow/srink 

E = 10     #Young modulus (everything is unit less since the simulation is not very physical anyways) Lower young modulus means more esay to deform means higher effect of the concentration
dE = 5     #difference in the young modulus as defined in the paper
c0 = 0.5   #reference concentration as defined in the paper
a0 = 1     #Inital radius of the tubes 
Nk = 20    #Quantity of softening agent added each step


#All those parameter interect in the following equation (a-a0)*E = alpha, for some constant.
#This equation being valid is equivalent to saying the elastic forces are constant meaning it has the time at each step to balance with the pressure force (this is of course a simplification)
#The pressure forces are not resolved because the hydrodynamics is far too complicared instead we just impose that the volume should be conserved
#In practice we try to minimize sum_i ((a_i-a0)-alpha/E_i)^2 under the constraint that the volume is conserved : sum_i a_i^2 = V

def RadiusEresolve(a,E): #This function takes as input the list of radius and the list of young modulus for each tube and return the solution of the minimization problem a list of radius na
    V = sum(a*a)         #The total volume
    h = np.sqrt(np.sum(np.square(np.full(len(E),a0)+alpha/E))/V)-1    #h is the lagrange multiplier imposing the conserved volume
    na = (np.full(len(E),a0)+alpha/E)/(1+h) 
    return na


# The other problem is how to represent the advection of the softening agent, in reality this is done thought cytoplasmic flows, but this is far outside the scope of this simulation 
# In pratice the effect is that the the fluid inside is advected in all direction pretty much equally the net effect of those is a kind of diffusion (This is NOT actual diffusion, one of the point of the paper is prescily to show that the relevant speed is the cytoplasmic speed and not the diffusion speed)
# The advection is done as self averaging effect, for each edge all the softening adgent is divided between the neigbouring edge (including itself) each contribution is weighted with respect to the radius squared since the flow rate is proportionnal to the area (the radius squared)

def advect(G,Gn,Ga):         #This function does what's described above for each edge we add to neighbouring edges a fraction of the number of softening agent.       
    n,m = len(G)//2,len(G[0])
    nGn = CreateEmptyGraph(m,n)
    for i in range(len(Gn)):
        for j in range(len(Gn[i])):
            if G[i][j]:
                P1,P2 = PointFromGraphIndice((i,j))
                Neigh = NeighbourVertice(G,P1,P2)
                neighV = np.sum([Ga[u][v]**2 for (u,v) in Neigh])
                for (u,v) in Neigh:
                    nGn[u][v] = nGn[u][v] + (Gn[i][j]*Ga[u][v]**2)/neighV
    return nGn

def updateRadius(G,Gn,Ga):        #This function takes Gn the number of solftening agent and the radius of each tube 
    n,m = len(G)//2,len(G[0])
    nGa = CreateEmptyGraph(m,n)
    Le = []
    La = []
    for i in range(len(Gn)):                 
        for j in range(len(Gn[i])):
            if G[i][j]:
                c = Gn[i][j]/Ga[i][j]**2           #we calculate the concentration as n/a^2
                Le.append(E-dE*c/(c+c0))           #This is the value of the actual young modulus using the formula from the paper
                La.append(Ga[i][j])
    nLa = RadiusEresolve(np.array(La),np.array(Le))  #Using the list of young modulus we update the radius of the tubes as described above
    k = 0
    for i in range(len(Gn)):                 #Put the new values in
        for j in range(len(Gn[i])):
            if G[i][j]:
                nGa[i][j] = nLa[k]
                k+=1
    return nGa

def Simulate(Nstep,G,Ga,Gn):                    #Now we put everyting together
    n,m = len(G)//2,len(G[0])
    L = [Ga]
    for _ in range(Nstep):                      #Nstep is the number of simulation step
        Gn[1][m//2] = Gn[1][m//2] + Nk          #we add at the origion a certain quantity of softening agent, this is the source of nutrient (could easly be modified to add different amount at different times)
        Gn = advect(G,Gn,Ga)                    #we advect the solftening agent
        Ga = updateRadius(G,Gn,Ga)              #we update the radius and then repeat 
        L.append(Ga)
    return L                                    #In the end we have a list of radius of each tube at each step

#The next functions are used for visualisation

def ColorGrade(t):     #this takes in a parameter t between -1 and 1 and transform's it into a color here we go from blue to gray to yellow (trying to mimick the color grading of the paper)
    g = 0.55
    if t < 0:
        return (g*(t+1),g*(t+1),-t+g*(t+1))
    else:
        return (g*(1-t)+t,g*(1-t)+0.8*t,g*(1-t))

def PatchCollection(G,Ga,cmult = 3):  # This function return a collection of patches (shapes) and colors that then will be drawn to the screen, cmult is a color multiplier, if one increase it's value the color are more intense
    n,m = len(G)//2,len(G[0])
    AvgA = np.sum([np.sum(Ga[i]*G[i]) for i in range(len(Ga))])/np.sum([np.sum(G[i]) for i in range(len(Ga))])  #We chose a color based on the distance to average radius
    patchs = []
    colors = []
    for i in range(len(Ga)):   #first the edges
        for j in range(len(Ga[i])):
            if G[i][j]:
                if i%2 == 0:
                    if j != m-1:
                        rect = ptc.Rectangle((5+13*j,4.5+13*i//2-Ga[i][j]/2),12,Ga[i][j])
                        colors.append(ColorGrade(np.tanh(cmult*(Ga[i][j]-AvgA))))  #the tanh function is used to clamp the value to -1,1
                        patchs.append(rect)
                else:
                    if i != 2*n-1:
                        rect = ptc.Rectangle((4.5+13*j-Ga[i][j]/2,13*i//2-1),Ga[i][j],12)
                        colors.append(ColorGrade(np.tanh(cmult*(Ga[i][j]-AvgA))))
                        patchs.append(rect)
    for i in range(m):       # then the vertices
        for j in range(n):
            neigh = neighbour(G,i,j)
            if len(neigh) > 1:        #we only plot a vertice if it has an edge connected to it
                a = []
                for (u,v) in neigh:
                    if (u,v) != (i,j):
                        ind = VerticeIndice(G,np.array([i,j]),np.array([u,v]))
                        a.append(np.tanh(cmult*(Ga[ind[0]][ind[1]]-AvgA)))       #the color of a vertice is the average color of the surounding edges
                t = np.mean(a)
                rect = ptc.Circle((4.5+13*i,4.5+13*j),1.5)
                colors.append(ColorGrade(t))
                patchs.append(rect)
    return collec.PatchCollection(patchs,animated=True),colors

def RenderAnimation(L,G,cmult = 3):   #This is the acutal function used to animate
    n,m = len(G)//2,len(G[0])
    fig, ax = plt.subplots()

    ax.set_axis_off()
    ax.set_aspect('equal', adjustable='box')

    def animate(frame):
        patchs,colors = PatchCollection(G,L[frame])
        ax.add_collection(patchs)
        patchs.set_facecolor(colors)
        return patchs,

    plt.ylim(0,13*n+3)
    plt.xlim(0,13*m+3)
    fig.tight_layout()
    ani = anim.FuncAnimation(fig,animate,frames = range(len(L)),interval = 300)  #interval is the time in ms between each image
    ani.save(pathlib.Path(path+"/animation.mp4"), dpi = 200)   #the animation is rendered to the file that the python file is in but it can be easly modified if needed (be careful since the all have the same name is there is a video file you want to save rename it or move it to a different file or it'll be deleted when a new one is created)

m,n = 9,18
G = GenerateRandomGraph(m,n,0.75)
G = ClearUnconnectedBits(G)
Gai = GenerateIniGa(G,a0)
Gn = CreateEmptyGraph(m,n)
L = Simulate(50,G,Gai,Gn)
RenderAnimation(L,G)

# decay part can simply be added by multiplying the n by a number less then 1 at each step


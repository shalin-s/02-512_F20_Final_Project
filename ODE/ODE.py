import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

## First ODEs

def ODE_1(k1, k2, k3, S0, I0, R0, D0, T, dt):
    inputs = (k1, k2, k3, dt)

    S = S0
    I = I0
    R = R0
    D = D0

    def simS(S, I, R, D, inputs, deltat, Scurr):
        (k1, k2, k3, dt) = inputs
        new = (Scurr + (deltat)*(-k1*S))
        return new

    def simI(S, I, R, D, inputs, deltat, Icurr):
        (k1, k2, k3, dt) = inputs
        new = (Icurr + (deltat)*(k1*S - (k2+k3)*I))
        return new

    def simR(S, I, R, D, inputs, deltat, Rcurr):
        (k1, k2, k3, dt) = inputs
        new = (Rcurr + (deltat)*(k2*I))
        return new

    def simD(S, I, R, D, inputs, deltat, Dcurr):
        (k1, k2, k3, dt) = inputs
        new = (Dcurr + (deltat)*(k3*I))
        return new

    yaxisS, yaxisI, yaxisR, yaxisD = [], [], [], []
    xaxis = []

    for t in range(T):
        xaxis.append(t)

        Shat = simS(S, I, R, D, inputs, dt/2, S)
        Ihat = simI(S, I, R, D, inputs, dt/2, I)
        Rhat = simR(S, I, R, D, inputs, dt/2, R)
        Dhat = simD(S, I, R, D, inputs, dt/2, D)

        newS = simS(Shat, Ihat, Rhat, Dhat, inputs, dt, S)
        newI = simI(Shat, Ihat, Rhat, Dhat, inputs, dt, I)
        newR = simR(Shat, Ihat, Rhat, Dhat, inputs, dt, R)
        newD = simD(Shat, Ihat, Rhat, Dhat, inputs, dt, D)

        S, I, R, D = newS, newI, newR, newD

        yaxisS.append(S)
        yaxisI.append(I)
        yaxisR.append(R)
        yaxisD.append(D)
    
    plt.plot(xaxis, yaxisS, color='skyblue')
    plt.plot(xaxis, yaxisI, color='olive')
    plt.plot(xaxis, yaxisR, color='palevioletred')
    plt.plot(xaxis, yaxisD, color='darkorange')
    plt.title(f"Fig. 1: S, I, R, and D over 1000 simulations")
    plt.xlabel("Time")
    plt.ylabel("People")

    Spatch = mpatches.Patch(color='skyblue', label=f'S (Final: {int(S)})')
    Ipatch = mpatches.Patch(color='olive', label=f'I (Final: {int(I)})')
    Rpatch = mpatches.Patch(color='palevioletred', label=f'R (Final: {int(R)})')
    Dpatch = mpatches.Patch(color='darkorange', label=f'D (Final: {int(D)})')
    plt.legend(handles=[Spatch,Ipatch,Rpatch,Dpatch])
    plt.show()

def run_1():
    k1 = 0.5
    k2 = 0.95
    k3 = 0.05
    total = 100
    I0 = 5
    S0 = total - I0
    R0 = 0
    D0 = 0
    T = 400
    dt = 0.1
    ODE_1(k1, k2, k3, S0, I0, R0, D0, T, dt)

def run_1_vacc():
    k1 = 0.5
    k2 = 0.95
    k3 = 0.05
    total = 100
    I0 = 5
    R0 = 30
    S0 = total - I0 - R0
    D0 = 0
    T = 400
    dt = 0.1
    ODE_1(k1, k2, k3, S0, I0, R0, D0, T, dt)

## Second ODEs

def ODE_2(k1, k2, k3, S0, I0, R0, D0, T, dt):
    inputs = (k1, k2, k3, dt)

    S = S0
    I = I0
    R = R0
    D = D0

    def simS(S, I, R, D, inputs, deltat, Scurr):
        (k1, k2, k3, dt) = inputs
        new = (Scurr + (deltat)*(-k1*S*I)/(S + I + R))
        return new

    def simI(S, I, R, D, inputs, deltat, Icurr):
        (k1, k2, k3, dt) = inputs
        new = (Icurr + (deltat)*((k1*S*I)/(S + I + R) - (k2+k3)*I))
        return new

    def simR(S, I, R, D, inputs, deltat, Rcurr):
        (k1, k2, k3, dt) = inputs
        new = (Rcurr + (deltat)*(k2*I))
        return new

    def simD(S, I, R, D, inputs, deltat, Dcurr):
        (k1, k2, k3, dt) = inputs
        new = (Dcurr + (deltat)*(k3*I))
        return new

    yaxisS, yaxisI, yaxisR, yaxisD = [], [], [], []
    xaxis = []

    for t in range(T):
        xaxis.append(t)

        Shat = simS(S, I, R, D, inputs, dt/2, S)
        Ihat = simI(S, I, R, D, inputs, dt/2, I)
        Rhat = simR(S, I, R, D, inputs, dt/2, R)
        Dhat = simD(S, I, R, D, inputs, dt/2, D)

        newS = simS(Shat, Ihat, Rhat, Dhat, inputs, dt, S)
        newI = simI(Shat, Ihat, Rhat, Dhat, inputs, dt, I)
        newR = simR(Shat, Ihat, Rhat, Dhat, inputs, dt, R)
        newD = simD(Shat, Ihat, Rhat, Dhat, inputs, dt, D)

        S, I, R, D = newS, newI, newR, newD

        yaxisS.append(S)
        yaxisI.append(I)
        yaxisR.append(R)
        yaxisD.append(D)
    
    plt.plot(xaxis, yaxisS, color='skyblue')
    plt.plot(xaxis, yaxisI, color='olive')
    plt.plot(xaxis, yaxisR, color='palevioletred')
    plt.plot(xaxis, yaxisD, color='darkorange')
    plt.title(f"Fig. 1: S, I, R, and D over 1000 simulations")
    plt.xlabel("Time")
    plt.ylabel("People")

    Spatch = mpatches.Patch(color='skyblue', label=f'S (Final: {int(S)})')
    Ipatch = mpatches.Patch(color='olive', label=f'I (Final: {int(I)})')
    Rpatch = mpatches.Patch(color='palevioletred', label=f'R (Final: {int(R)})')
    Dpatch = mpatches.Patch(color='darkorange', label=f'D (Final: {int(D)})')
    plt.legend(handles=[Spatch,Ipatch,Rpatch,Dpatch])
    plt.show()

def run_2():
    k1 = 0.5
    k2 = 0.95
    k3 = 0.05
    total = 100
    I0 = 5
    S0 = total - I0
    R0 = 0
    D0 = 0
    T = 400
    dt = 0.1
    ODE_2(k1, k2, k3, S0, I0, R0, D0, T, dt)

def run_2_vacc():
    k1 = 0.5
    k2 = 0.95
    k3 = 0.05
    total = 100
    I0 = 5
    R0 = 30
    S0 = total - I0 - R0
    D0 = 0
    T = 400
    dt = 0.1
    ODE_2(k1, k2, k3, S0, I0, R0, D0, T, dt)

## Third ODE (like first but with masks)

def ODE_3(k1a, k1b, k2, k3, Sa0, Sb0, I0, R0, D0, T, dt):
    inputs = (k1a, k1b, k2, k3, dt)

    Sa = Sa0
    Sb = Sb0
    I = I0
    R = R0
    D = D0

    def simSa(Sa, Sb, I, R, D, inputs, deltat, Sacurr):
        (k1a, k1b, k2, k3, dt) = inputs
        new = (Sacurr + (deltat)*(-k1a*Sa))
        return new
    
    def simSb(Sa, Sb, I, R, D, inputs, deltat, Sbcurr):
        (k1a, k1b, k2, k3, dt) = inputs
        new = (Sbcurr + (deltat)*(-k1b*Sb))
        return new

    def simI(Sa, Sb, I, R, D, inputs, deltat, Icurr):
        (k1a, k1b, k2, k3, dt) = inputs
        new = (Icurr + (deltat)*(k1a*Sa + k1b*Sb - (k2+k3)*I))
        return new

    def simR(Sa, Sb, I, R, D, inputs, deltat, Rcurr):
        (k1a, k1b, k2, k3, dt) = inputs
        new = (Rcurr + (deltat)*(k2*I))
        return new

    def simD(Sa, Sb, I, R, D, inputs, deltat, Dcurr):
        (k1a, k1b, k2, k3, dt) = inputs
        new = (Dcurr + (deltat)*(k3*I))
        return new

    yaxisSa, yaxisSb, yaxisI, yaxisR, yaxisD = [], [], [], [], []
    xaxis = []

    for t in range(T):
        xaxis.append(t)

        Sahat = simSa(Sa, Sb, I, R, D, inputs, dt/2, Sa)
        Sbhat = simSb(Sa, Sb, I, R, D, inputs, dt/2, Sb)
        Ihat = simI(Sa, Sb, I, R, D, inputs, dt/2, I)
        Rhat = simR(Sa, Sb, I, R, D, inputs, dt/2, R)
        Dhat = simD(Sa, Sb, I, R, D, inputs, dt/2, D)

        newSa = simSa(Sahat, Sbhat, Ihat, Rhat, Dhat, inputs, dt, Sa)
        newSb = simSb(Sahat, Sbhat, Ihat, Rhat, Dhat, inputs, dt, Sb)
        newI = simI(Sahat, Sbhat, Ihat, Rhat, Dhat, inputs, dt, I)
        newR = simR(Sahat, Sbhat, Ihat, Rhat, Dhat, inputs, dt, R)
        newD = simD(Sahat, Sbhat, Ihat, Rhat, Dhat, inputs, dt, D)

        Sa, Sb, I, R, D = newSa, newSb, newI, newR, newD

        yaxisSa.append(Sa)
        yaxisSb.append(Sb)
        yaxisI.append(I)
        yaxisR.append(R)
        yaxisD.append(D)
    
    plt.plot(xaxis, yaxisSa, color='skyblue')
    plt.plot(xaxis, yaxisSb, color='royalblue')
    plt.plot(xaxis, yaxisI, color='olive')
    plt.plot(xaxis, yaxisR, color='palevioletred')
    plt.plot(xaxis, yaxisD, color='darkorange')
    plt.title(f"Fig. 2: Sa, Sb, I, R, and D over 1000 simulations (30% vacc.)")
    plt.xlabel("Time")
    plt.ylabel("People")

    Sapatch = mpatches.Patch(color='skyblue', label=f'Sa (Final: {int(Sa)})')
    Sbpatch = mpatches.Patch(color='royalblue', label=f'Sb (Final: {int(Sb)})')
    Ipatch = mpatches.Patch(color='olive', label=f'I (Final: {int(I)})')
    Rpatch = mpatches.Patch(color='palevioletred', label=f'R (Final: {int(R)})')
    Dpatch = mpatches.Patch(color='darkorange', label=f'D (Final: {int(D)})')
    plt.legend(handles=[Sapatch,Sbpatch,Ipatch,Rpatch,Dpatch])
    plt.show()

def run_3():
    k1a = 0.1
    k1b = 0.5
    k2 = 0.95
    k3 = 0.05
    total = 100
    I0 = 5
    S0 = total - I0
    Sa0 = 0.3*S0
    Sb0 = 0.7*S0
    R0 = 0
    D0 = 0
    T = 400
    dt = 0.1
    ODE_3(k1a, k1b, k2, k3, Sa0, Sb0, I0, R0, D0, T, dt)

def run_3_vacc():
    k1a = 0.1
    k1b = 0.5
    k2 = 0.95
    k3 = 0.05
    total = 100
    I0 = 5
    R0 = 30
    S0 = total - I0 - R0
    Sa0 = 0.3*S0
    Sb0 = 0.7*S0
    D0 = 0
    T = 400
    dt = 0.1
    ODE_3(k1a, k1b, k2, k3, Sa0, Sb0, I0, R0, D0, T, dt)
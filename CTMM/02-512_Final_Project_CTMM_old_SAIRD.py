# Code for 02-512 Final Project, Fall 2020, by Shalin Shah.
# Tested in Python 3.6.5

import sys
import numpy as np
import numpy.random
from datetime import datetime
from matplotlib import pyplot as plt
np.seterr(divide='ignore')

def log_print(text):
    print("{}: {}".format(datetime.now(), text))

# This article is helpful for standard terminology on these sorts of models: https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology

# indices in the state array
Q_S = 0
Q_A = 1
Q_I = 2
Q_R = 3
Q_D = 4
Q_N = 5
NUM_STATE_VARS = 6

STATE_NAMES = [None] * NUM_STATE_VARS
STATE_NAMES[Q_S] = "Susceptible"
STATE_NAMES[Q_A] = "Asymptomatic"
STATE_NAMES[Q_I] = "Infectious"
STATE_NAMES[Q_R] = "Recovered"
STATE_NAMES[Q_D] = "Dead"
STATE_NAMES[Q_N] = "Total #People Alive"

STATE_COLORS = [None] * NUM_STATE_VARS
STATE_COLORS[Q_S] = "blue"
STATE_COLORS[Q_A] = "yellow"
STATE_COLORS[Q_I] = "red"
STATE_COLORS[Q_R] = "green"
STATE_COLORS[Q_D] = "purple"
STATE_COLORS[Q_N] = "black"

# indices in the array of possible transitions at each step
T_ASYMPTOMATIC = 0
T_INFECTIOUS = 1
T_RECOVER_ASYMPTOMATIC = 2
T_RECOVER_INFECTIOUS = 3
T_DEATH = 4
NUM_EVENT_TYPES = 5

# Hardcoded parameter values (Ls are rates/lambdas)
N0 = 10000 # Initial Number of people
A0 = 10 # Initial Number of asymptomatic infected people
L_ASYMPTOMATIC_SPREAD = 0.12
L_INFECTIOUS_SPREAD = 0.09
L_INCUBATION = 0.1
L_RECOVER_ASYMPTOMATIC = 0.05
L_RECOVER_INFECTIOUS = 0.09
L_DEATH = 0.01


# Uncomment this block to read the parameter values from a file instead of them being hardcoded
'''
INPUT_FILE = sys.argv[1]
original_lines = open(INPUT_FILE).read().splitlines()
N0 = int(original_lines[0])
A0 = int(original_lines[1])
L_ASYMPTOMATIC_SPREAD = float(original_lines[3])
L_INFECTIOUS_SPREAD = float(original_lines[4])
L_INCUBATION = float(original_lines[5])
L_RECOVER = float(original_lines[6])
L_DEATH = float(original_lines[7])
'''

t_curr = 0

Q = np.zeros(NUM_STATE_VARS).astype(np.int64)
Q[Q_N] = N0
Q[Q_S] = N0 - A0
Q[Q_A] = A0
Q[Q_I] = 0 # redundant, but here for clarity.
Q[Q_R] = 0 # redundant, but here for clarity.

T = np.zeros(NUM_EVENT_TYPES)

QS_list = []
QS_list.append(Q.copy())

times_list = []
times_list.append(t_curr)

iter_count = 0
while (max(Q[Q_A], Q[Q_I]) > 0):
	if (iter_count % 100 == 0):
		log_print("Iteration {}, Time {}, State ([S, A, I, R, D, N]) {}".format(iter_count, t_curr, Q))

	lambdas = np.zeros(T.shape)
	lambdas[T_ASYMPTOMATIC] = (L_ASYMPTOMATIC_SPREAD * Q[Q_A] * (Q[Q_S] / Q[Q_N])) + (L_INFECTIOUS_SPREAD * Q[Q_I] * (Q[Q_S] / Q[Q_N])) # New infection (always asymptomatic at first)
	lambdas[T_RECOVER_ASYMPTOMATIC] = L_RECOVER_ASYMPTOMATIC * Q[Q_A] # Recover from asymptomatic infection
	lambdas[T_INFECTIOUS] = L_INCUBATION * Q[Q_A] # Progression of infection from asymptomatic to symptomatic ("infectious")
	lambdas[T_RECOVER_INFECTIOUS] = L_RECOVER_INFECTIOUS * Q[Q_I] # Recover from asymptomatic infection
	lambdas[T_DEATH] = L_DEATH * Q[Q_I] # Die from (necessarily symptomatic) infection
	
	betas = 1/lambdas # inf handling (1/0 cases) is handled properly automatically

	T = np.random.exponential(scale=betas)

	min_time_idx = np.argmin(T)

	t_curr += T[min_time_idx] # Update time no matter what

	if (min_time_idx == T_ASYMPTOMATIC):
		Q[Q_S] -= 1
		Q[Q_A] += 1

	elif (min_time_idx == T_INFECTIOUS):
		Q[Q_A] -= 1
		Q[Q_I] += 1

	elif (min_time_idx == T_RECOVER_ASYMPTOMATIC):
		Q[Q_A] -= 1
		Q[Q_R] += 1

	elif (min_time_idx == T_RECOVER_INFECTIOUS):
		Q[Q_I] -= 1
		Q[Q_R] += 1

	elif (min_time_idx == T_DEATH):
		Q[Q_I] -= 1
		Q[Q_D] += 1
		Q[Q_N] -= 1

	QS_list.append(Q.copy())
	times_list.append(t_curr)

	iter_count += 1


QS = np.array(QS_list).T # each row is the values for that state

# Print final results and other relevant information
log_print("")
log_print("Final State ([S, A, I, R, D, N]): {}".format(Q))

for Q_idx in range(NUM_STATE_VARS):
	plt.plot(times_list, QS[Q_idx, :], color=STATE_COLORS[Q_idx], label=STATE_NAMES[Q_idx] + " (Final: {})".format(Q[Q_idx]))
plt.legend()
plt.title("CTMM Simulation of SAIRD Model ({} initial asymptomatic infections, among {} people)".format(A0, N0))
plt.xlabel("Time")
plt.ylabel("Number of People")
fig = plt.gcf()
fig.set_size_inches(9, 9, forward=True)
plt.savefig("CTMM_SAIRD_1.png")
plt.close()

log_print("Done")


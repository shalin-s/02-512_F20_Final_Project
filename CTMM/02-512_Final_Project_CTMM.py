# Code for 02-512 Final Project, Fall 2020, by Shalin Shah.
# Tested in Python 3.6.5

# WORK IN PROGRESS, THIS IS NOT FULLY FINISHED YET

import sys
import numpy as np
import numpy.random
from datetime import datetime
from matplotlib import pyplot as plt
np.seterr(divide='ignore')

def log_print(text):
    print("{}: {}".format(datetime.now(), text))

# This article is helpful for standard terminology on these sorts of models: https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology

# rows in the state array
Q_S = 0
Q_I = 1
Q_R = 2
Q_D = 3
Q_N = 4
NUM_STATE_VARS = 5

QC_BASE = 0
QC_MASK = 1
NUM_STATE_COLS = 2
#WAIT, SHOULD I REPRESENT Q AS 5X2, OR 9X1?? WHICH WOULD BE EASIER TO DEAL WITH? I feel like 5x2 would be better logically, but might be kind of annoying to add columns everywhere
#OR should Q be 2x5? Idk yet.

STATE_NAMES = [None] * NUM_STATE_VARS
STATE_NAMES[Q_S] = "Susceptible"
STATE_NAMES[Q_I] = "Infectious"
STATE_NAMES[Q_R] = "Recovered"
STATE_NAMES[Q_D] = "Dead"
STATE_NAMES[Q_N] = "Total #People Alive"

STATE_COLORS = [None] * NUM_STATE_VARS
STATE_COLORS[Q_S] = "blue"
STATE_COLORS[Q_I] = "red"
STATE_COLORS[Q_R] = "green"
STATE_COLORS[Q_D] = "purple"
STATE_COLORS[Q_N] = "black"

# indices in the array of possible transitions at each step
T_INFECTION = 0
T_RECOVER = 1
T_DEATH = 2
NUM_EVENT_TYPES = 3

# Hardcoded parameter values (Ls are rates/lambdas)
N0 = 10000 # Initial Number of people
I0 = 10 # Initial Number of infected people
L_INFECTION = 2.5
L_RECOVER = 0.95
L_DEATH = 0.05


# Uncomment this block to read the parameter values from a file instead of them being hardcoded
'''
INPUT_FILE = sys.argv[1]
original_lines = open(INPUT_FILE).read().splitlines()
N0 = int(original_lines[0])
I0 = int(original_lines[1])
L_INFECTIOUS_SPREAD = float(original_lines[2])
L_RECOVER = float(original_lines[3])
L_DEATH = float(original_lines[4])
'''

t_curr = 0

Q = np.zeros(NUM_STATE_VARS).astype(np.int64)
Q[Q_N] = N0
Q[Q_S] = N0 - I0
Q[Q_I] = I0
Q[Q_R] = 0 # redundant, but here for clarity.
Q[Q_D] = 0 # redundant, but here for clarity.

T = np.zeros(NUM_EVENT_TYPES)

QS_list = []
QS_list.append(Q.copy())

times_list = []
times_list.append(t_curr)

iter_count = 0
while (Q[Q_I] > 0):
	if (iter_count % 100 == 0):
		log_print("Iteration {}, Time {}, State ([S, I, R, D, N]) {}".format(iter_count, t_curr, Q))

	lambdas = np.zeros(T.shape)
	lambdas[T_INFECTION] = L_INFECTION * Q[Q_I] * (Q[Q_S] / Q[Q_N]) # New infection
	lambdas[T_RECOVER] = L_RECOVER * Q[Q_I] # Recover from infection
	lambdas[T_DEATH] = L_DEATH * Q[Q_I] # Die from infection
	
	betas = 1/lambdas # inf handling (1/0 cases) is handled properly automatically

	T = np.random.exponential(scale=betas)

	min_time_idx = np.argmin(T)

	t_curr += T[min_time_idx] # Update time no matter what

	if (min_time_idx == T_INFECTION):
		Q[Q_S] -= 1
		Q[Q_I] += 1

	elif (min_time_idx == T_RECOVER):
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
log_print("Final State ([S, I, R, D, N]): {}".format(Q))

for Q_idx in range(NUM_STATE_VARS):
	plt.plot(times_list, QS[Q_idx, :], color=STATE_COLORS[Q_idx], label=STATE_NAMES[Q_idx] + " (Final: {})".format(Q[Q_idx]))
plt.legend()
plt.title("CTMM Simulation of SIRD Model ({} initial asymptomatic infections, among {} people)".format(I0, N0))
plt.xlabel("Time")
plt.ylabel("Number of People")
fig = plt.gcf()
fig.set_size_inches(9, 9, forward=True)
plt.savefig("CTMM_SIRD_Masks_1.png")
plt.close()

log_print("Done")


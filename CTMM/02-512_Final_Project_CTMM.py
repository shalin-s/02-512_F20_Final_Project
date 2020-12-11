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

# This article is helpful for reference and standard terminology on SIR and similar models, including SIRD): https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology

# rows in the state array
Q_S = 0
Q_I = 1
Q_R = 2
Q_D = 3
Q_N = 4
NUM_STATE_VARS = 5

# Columns in the state array (mask and base (no masks))
QC_BASE = 0
QC_MASK = 1
NUM_STATE_COLS = 2

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

STATE_COL_NAMES = [None] * NUM_STATE_COLS
STATE_COL_NAMES[QC_BASE] = "No Mask"
STATE_COL_NAMES[QC_MASK] = "Mask"

STATE_COL_LINESTYLES = [None] * NUM_STATE_COLS
STATE_COL_LINESTYLES[QC_BASE] = "-"
STATE_COL_LINESTYLES[QC_MASK] = "--"

# indices in the array of possible transition types at each step (events)
T_INFECTION = 0
T_RECOVER = 1
T_DEATH = 2
NUM_EVENT_TYPES = 3

# Hardcoded parameter values (L-values are rates, i.e. lambdas):
N0 = 10000 # Initial Number of people
MASK_PROPORTION = 0.3 # Please make sure that MASK_PROPORTION * I0 is an integer, for consistency
I0 = 10 # Initial Number of infected people
LS_INFECTION = np.zeros((NUM_STATE_COLS, NUM_STATE_COLS))
LS_INFECTION[QC_BASE, QC_BASE] = 5.0 # infection spread rate from base to base
LS_INFECTION[QC_BASE, QC_MASK] = 1.0 # infection spread rate from base to mask
LS_INFECTION[QC_MASK, QC_BASE] = 1.0 # infection spread rate from mask to base
LS_INFECTION[QC_MASK, QC_MASK] = 0.5 # infection spread rate from mask to mask
L_RECOVER = 0.95
L_DEATH = 0.05


# Uncomment this block to read the parameter values from a file instead of them being hardcoded
'''
INPUT_FILE = sys.argv[1]
original_lines = open(INPUT_FILE).read().splitlines()
N0 = int(original_lines[0])
I0 = int(original_lines[1])
MASK_PROPORTION = float(original_lines[2])
LS_INFECTION[QC_BASE, QC_BASE] = float(original_lines[3])
LS_INFECTION[QC_BASE, QC_MASK] = float(original_lines[4])
LS_INFECTION[QC_MASK, QC_BASE] = float(original_lines[5])
LS_INFECTION[QC_MASK, QC_MASK] = float(original_lines[6])
L_RECOVER = float(original_lines[7])
L_DEATH = float(original_lines[8])
'''

t_curr = 0

Q = np.zeros((NUM_STATE_VARS, NUM_STATE_COLS)).astype(np.int64)
Q[Q_N, QC_MASK] = int(round(N0 * MASK_PROPORTION))
Q[Q_N, QC_BASE] = N0 - Q[Q_N, QC_MASK]
Q[Q_S, QC_MASK] = int(round((N0 - I0) * MASK_PROPORTION))
Q[Q_S, QC_BASE] = (N0 - I0) - Q[Q_S, QC_MASK]
Q[Q_I, QC_MASK] = int(round((I0) * MASK_PROPORTION))
Q[Q_I, QC_BASE] = I0 - Q[Q_I, QC_MASK]
Q[Q_R, :] = 0 # redundant, but here for clarity.
Q[Q_D, :] = 0 # redundant, but here for clarity.


# This is only for testing of vaccination scenario (move some people to "recovered" state initially)
'''
vaccination_rate = 0.3
num_vaccinated = vaccination_rate * N0

Q[Q_S, QC_MASK] = int(round((N0 - I0 - num_vaccinated) * MASK_PROPORTION))
Q[Q_S, QC_BASE] = (N0 - I0 - num_vaccinated) - Q[Q_S, QC_MASK]

Q[Q_R, QC_MASK] = int(round(num_vaccinated * MASK_PROPORTION))
Q[Q_R, QC_BASE] = num_vaccinated - Q[Q_R, QC_MASK]

assert(np.sum(Q[Q_N]) == N0)
'''

T = np.zeros((NUM_EVENT_TYPES, NUM_STATE_COLS))

QS_list = [] # state history
QS_list.append(Q.copy())

times_list = []
times_list.append(t_curr)

iter_count = 0
while (np.sum(Q[Q_I, :]) > 0):
	if (iter_count % 100 == 0):
		log_print("Iteration {}, Time {}, No Mask State ([S, I, R, D, N]): {}".format(iter_count, t_curr, Q[:, QC_BASE]))
		log_print("Iteration {}, Time {}, Mask State ([S, I, R, D, N]): {}".format(iter_count, t_curr, Q[:, QC_MASK]))

	lambdas = np.zeros(T.shape)
	lambdas[T_INFECTION, QC_BASE] = LS_INFECTION[QC_MASK, QC_BASE] * Q[Q_I, QC_MASK] * (Q[Q_S, QC_BASE] / np.sum(Q[Q_N])) + LS_INFECTION[QC_BASE, QC_BASE] * Q[Q_I, QC_BASE] * (Q[Q_S, QC_BASE] / np.sum(Q[Q_N]))
	lambdas[T_INFECTION, QC_MASK] = LS_INFECTION[QC_MASK, QC_MASK] * Q[Q_I, QC_MASK] * (Q[Q_S, QC_MASK] / np.sum(Q[Q_N])) + LS_INFECTION[QC_BASE, QC_MASK] * Q[Q_I, QC_BASE] * (Q[Q_S, QC_MASK] / np.sum(Q[Q_N]))
	
	lambdas[T_RECOVER] = L_RECOVER * Q[Q_I] # Recover from infection (this works elementwise)
	lambdas[T_DEATH] = L_DEATH * Q[Q_I] # Die from infection (this works elementwise)
	
	betas = 1/lambdas # inf handling (1/0 cases) is handled properly automatically

	T = np.random.exponential(scale=betas)

	min_i, min_j = tuple(np.unravel_index(np.argmin(T, axis=None), T.shape))

	t_curr += T[min_i, min_j] # Update time no matter what

	if (min_i == T_INFECTION):
		Q[Q_S, min_j] -= 1
		Q[Q_I, min_j] += 1

	elif (min_i == T_RECOVER):
		Q[Q_I, min_j] -= 1
		Q[Q_R, min_j] += 1

	elif (min_i == T_DEATH):
		Q[Q_I, min_j] -= 1
		Q[Q_D, min_j] += 1
		Q[Q_N, min_j] -= 1

	QS_list.append(Q.copy())
	times_list.append(t_curr)

	iter_count += 1

QS = np.array(QS_list)
print(QS.shape)

# Print final results
log_print("")
log_print("Final No Mask State ([S, I, R, D, N]): {}".format(Q[:, QC_BASE]))
log_print("Final Mask State ([S, I, R, D, N]): {}".format(Q[:, QC_MASK]))

# Everything from here on out is just building some nice-looking plots
for Q_col in range(NUM_STATE_COLS):
	for Q_state in range(NUM_STATE_VARS):
		label="{} ({}) (Final: {})".format(STATE_NAMES[Q_state], STATE_COL_NAMES[Q_col], Q[Q_state, Q_col])
		plt.plot(times_list, QS[:, Q_state, Q_col], color=STATE_COLORS[Q_state], linestyle=STATE_COL_LINESTYLES[Q_col], label=label)

figure_text = ""
figure_text += "Mask Death Rate: {}/{} = {:0.3f}%. ".format(Q[Q_D, QC_MASK], QS[0, Q_N, QC_MASK], 100 * Q[Q_D, QC_MASK]/QS[0, Q_N, QC_MASK]) 
figure_text += "No-Mask Death Rate: {}/{} = {:0.3f}%.\n".format(Q[Q_D, QC_BASE], QS[0, Q_N, QC_BASE], 100 * Q[Q_D, QC_BASE]/QS[0, Q_N, QC_BASE]) 
figure_text += "Overall Death Rate: {}/{} = {:0.3f}%. ".format(N0 - np.sum(Q[Q_N]), N0, 100 * (N0 - np.sum(Q[Q_N]))/ N0)
plt.figtext(0.5, 0.01, figure_text, wrap=True, horizontalalignment='center', fontsize=10)

plt.legend()
plt.title("CTMM Simulation of SIRD Model ({}/{} initial infections, {}% mask adoption)".format(I0, N0, MASK_PROPORTION*100))
plt.xlabel("Time")
plt.ylabel("Number of People")
fig = plt.gcf()
fig.set_size_inches(9, 9, forward=True)
plt.savefig("CTMM_SIRD_Masks_2.png")
plt.close()

log_print("Done")

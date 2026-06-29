from z3 import *

# this code was basically me learning how to learn Z3, it's not anything useful!

s = Solver()

LacI = Bool('LacI')
TetR = Bool('TetR')

s.add(Implies(LacI, Not(TetR)))

s.add(Implies(TetR, Not(LacI)))

s.add(LacI == True)
s.add(TetR == True)

result = s.check()

if result == sat:
    print("Circuit is logically consistent!")
    print(s.model())
else:
    print("Circuit is INCONSISTENT (Biological Paradox Found).")
    print("Unsat Core:", s.unsat_core())

# Define a function: Is the Promoter ON? 
# Input: Repressor Protein (Real), Threshold (Real)
# Output: Boolean
def RepressiblePromoter(protein_conc, threshold):
    return protein_conc < threshold

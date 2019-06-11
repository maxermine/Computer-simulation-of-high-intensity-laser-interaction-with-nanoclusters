from main import * # module main connection

def main(): # main function
 Q1 = 1 # first ball charge [relative units]
 Q2 = 1 # second ball charge [relative units]
 R1 = 1 # first ball radius [relative units]
 R2 = 1 # second ball radius [relative units]
 DX = 2 # distance between balls [relative units]
 N = 10 # number of nodes in radius, phi, tetta []

 f = forceComputeByParameters(Q1, Q2, R1, R2, DX, N) # compute the interaction force between two balls by parameters [relative units]

 print(DX, f) # output the distance and force

main()

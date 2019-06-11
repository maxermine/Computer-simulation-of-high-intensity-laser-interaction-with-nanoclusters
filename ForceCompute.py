from main import * # module main connection

ne = 5.5E+27 # electron cluster density [1 / m3]
ni = 5.5E+27 # ion cluster density [1 / m3]

Re = 30E-9 # electron cluster radius [m]
Ri = 30E-9 # ion cluster radius [m]

N = 10 # number of nodes in radius, phi, tetta []

def main(): # main function
 Qe = e * ne * 4 * math.pi * Re ** 3 / 3 # total electron cluster charge [C]
 Qi = e * ni * 4 * math.pi * Ri ** 3 / 3 # total ion cluster charge [C]

 x = 10E-9 # distance between clusters [m]

 f = forceComputeByParameters(Qe, Qi, Re, Ri, x, N) / (4 * math.pi * eo) # compute interaction force between clusters [N]

 print(x, f) # output the distance and force

main()

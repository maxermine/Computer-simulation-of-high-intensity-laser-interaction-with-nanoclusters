from main import * # module main connection

ne = 5.5E+27 # electron cluster density [1 / m3]
ni = 5.5E+27 # ion cluster density [1 / m3]

Re = 30E-9 # electron cluster radius [m]
Ri = 30E-9 # ion cluster radius [m]

def main(): # main function
 file = "F.dat" # table file name

 T = readFile(file) # read table from file to memory

 x = 10E-9 # distance between clusters [m]

 f = forceComputeByTable(ne, ni, Re, Ri, x, T) # compute interaction force between clusters using preculculated table [N]

 print(x, f) # output the distance and force

main()

from main import * # module main connection

ne = 27E+27 # electron cluster density [1 / m3]
ni = 27E+27 # ion cluster density [1 / m3]

Re = 1.39E-9 # electron cluster radius [m]
Ri = 1.39E-9 # ion cluster radius [m]

b = 1.83E+15 # oscillation decay parameter [1 / s]

def Laser(t): # laser function, parameters: t time [s]
 Qe = e * ne * 4 * math.pi * Re ** 3 / 3 # electron cluster total charge [C]

 I = 1E+16 # laser radiation intensity [W / m2]
 E = math.sqrt(2 * I / (c * eo)) # laser electric field strength [V / m]

 to = 50E-15 # laser pulse width [s]
 L = 4.36E-7 # laser wave length [m]
 w = 2 * math.pi * c / L # cyclic laser oscillation frequency [1 / s] 

 return E * Qe * math.cos(w * t) * math.cos(0.5 * math.pi * t ** 2 / to ** 2) # force interaction between electronic cluster and laser electric field [N]

def forceCompute(t, v, x, T): # electron cluster acceleration for Newton's differential equation, t time [s], v velocity [m / s], x coordinate [m], T table preculculated force
 Me = me * ne * 4 * math.pi * Re ** 3 / 3 # electron cluster total mass [kg]

 return (Laser(t) - forceComputeByTable(ne, ni, Re, Ri, x, T)) / Me - b * v # electron cluster acceleration in laser wave electric field [m / s2]

def main(): # main function
 file = "F.dat" # table file name

 T = readFile(file) # read table from file to memory

 sw = 4 # print switch: 1 electron cluster coordinate x / Re [], 2 electron cluster velocity [m / s], 3 electron cluster acceleration [m / s2], 4 electron cluster temperature [ev]

 t = -50E-15 # laser start time [s]
 x = 0 # start electron cluster coordinate [m]
 v = 0 # start electron cluster velocity [m / s]

 tm = 50E-15 # laser end time [s]
 h = tm / 1000 # iteration time step [s] 

 RungeKutta(tm, h, sw, forceCompute, t, v, x, T, Re, b) # Runge-Kutta solve method for differential equation

main()

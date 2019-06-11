import math # module math connection
from array import * # module array connection

e = 1.602176620898E-19 # electron charge [C]
me = 9.1093835611E-31 # electron mass [kg]
eo = 8.85418781762039E-12 # electrical constant [F / m]
c = 299792458 # velocity of light [m / s]

def forceComputeByParameters(Q1, Q2, R1, R2, DX, N): # compute interaction force between two balls for parameters: Q1 Q2 charges of balls, R1 R2 radius of balls, DX distance between balls, N number of nodes
 fx = 0 # force initial value

 dr = R1 / N # step for radius
 dfi = 2 * math.pi / N # step for fi
 dth = math.pi / N # step for tetta

 ro = Q1 / (4 * math.pi * R1 ** 3 / 3) # charge density of the first ball

 fi = dfi / 2 # initial value for fi angle in the middle of the cell

 while fi < (2 * math.pi): # cycle for angle fi
  th = dth / 2 # # initial value for tetta angle in the middle of the cell

  while th < math.pi: # cycle for angle tetta
   r = dr / 2 # # initial value for radius in the middle of the cell

   while r < R1: # cycle for radius
    x = r * math.sin(th) * math.cos(fi) # x coordinate in spherical coordinate system   
    y = r * math.sin(th) * math.sin(fi) # y coordinate in spherical coordinate system 
    z = r * math.cos(th) # z coordinate in spherical coordinate system 

    q = ro * r ** 2 * dr * math.sin(th) * dth * dfi # charge for small cell in the first ball

    rx = DX - x # distance along the x axis between the small cell in the first ball and the center of the second ball
    d = math.sqrt(rx ** 2 + y ** 2 + z ** 2) # full distance between the small cell in the first ball and the center of the second ball

    if d < R2: # if the distance is less than the radius of the second ball, then a small cell is inside the second ball
     fx += Q2 * q * rx / (R2 ** 3) # the interaction force is calculated by the Gauss theorem and summed
    else: # if the distance is greater than the radius of the second ball, then a small cell is outside the second ball
     fx += Q2 * q * rx / (d ** 3) # the interaction force is calculated using Coulomb formula and summed 

    r += dr # next cell for radius
   th += dth # next cell for tetta angle
  fi += dfi # next cell for fi angle

 return fx # the final value of the interaction force between the balls [relative units]

def forceComputeByTable(n1, n2, r1, r2, d, T): # compute the force of interaction between two balls using preculculated table for parameters: n1 n2 charge density of balls [1 / m3], r1 r2 radius of balls [m], d distance between the balls [m], T preculculated table
 Q1 = e * n1 * 4 * math.pi * r1 ** 3 / 3 # total first ball charge [C]
 Q2 = e * n2 * 4 * math.pi * r2 ** 3 / 3 # total second ball charge [C]

 ratio = r1 / r2 # ratio of balls radius 

 r = r2 # scale factor equated to second ball radius
 if r1 > r: # first ball is large than second ball
  r = r1 # scale factor equated to first ball radius
 
 x = d / r # the relative distance between the balls through a scale factor

 z = 1 # force projection factor is first taken as positive
 if x < 0: # if the one ball has shifted to the negative side of the coordinate axis so the coordinate is negative
  x = -x # translate coordinate to positive distance
  z = -1 # force projection factor is taken as negative

 return z * Q1 * Q2 * forceFromTable(x, ratio, T) / (4 * math.pi * eo * r ** 2) # force of the interaction between two charged balls [N]

def forceFromTable(x, ratio, t): # function receives the normalized force from the table or makes a linear approximation if there are no nodes in the table for parameters: x distance between balls, ratio radius ratio for balls, t table with preculculated forces
 if ratio > 1: # radius ratio is such that the second ball is smaller than the first
  ratio = 1. / ratio # swap balls to make the second ball larger than the first

 if x >= (1 + ratio): # the distance between the balls is such that the balls do not intersect
  return 1. / x ** 2 # return force through Coulomb's law formula

 a_ = t[0] # first table entry for radius ratio
 b_ = t[0] # first table entry for distance between the balls
 c_ = t[0] # first table entry for given distance between the balls
 d_ = t[0] # first table entry for given radius ratio

 for a in t: # table loop
  if a[0] == x and a[2] == ratio: # distance and radius ratio are matched
   return a[1] # return force value from table
 
  if a[2] == ratio: # ratio of the radius are coincided
   if a[0] > x: # distance is out of table
    return a_[1] + dF(a[1], a_[1], a[0], a_[0], x) # returns force value through linear approximation
   
   a_ = a # memorized entry with given radius ratio value

  if a[0] == x: # distance between balls is matched
   if a[2] < ratio: # radius ratio is out of table
    return a[1] + dF(b_[1], a[1], b_[2], a[2], ratio) # return force value through linear approximation

   b_ = a # memorized entry with given distance value

  if a[2] > ratio: # radius ratio is still large
   if a[0] < x: # the distance is still small
    d_ = a # memorized record with these values

  if a[2] < ratio: # ratio of radius is already less
   if a[0] > x: # distance is already large
    return c_[1] + dF(a[1], c_[1], a[0], c_[0], x) + dF(d_[1], c_[1], d_[2], c_[2], ratio) # return force value through linear approximation

  c_ = a # remember the current table entry

 return -1 # return error code if the algorithm did not work properly

def dF(y2, y1, x2, x1, x): # calculate increment for function through derivative of the function and the increment of the argument for parameters: y2 y1 function values, x2 x1 argument values, x point between x2 and x1
 if x2 == x1: # this is one point
  return 0 # return zero increment
 else: # this is two different points
  return (x - x1) * (y2 - y1) / (x2 - x1) # return function incriment 

def readFile(file): # read table from file into memory
 t = [] # initiate list for table entries
 with open(file) as f: # open file for reading
  for line in f: # read a record from file
   t.append([float(x) for x in line.split()]) # add fields from record to list

 return t # return address to the list

def RungeKutta(tm, h, sw, forceCompute, t, v, x, T, Re, b): # Runge-Kutta solve method for differential equation with parameters: tm end time, h step time, sw print switch, forceCompute user predefined function for culculation acceleration, t start time, v start velocity, x start coordinate, T preculculated table with forces, Re electron cluster radius, b oscilation decay parameter
 S = 0 # initial value for temperature integral

 while t < tm: # time cycle
  if sw == 1: # switch to print the relative coordinate
   print(t * 1E+15, x / Re) # output time [fs] and electron cluster coordinate [relative units]

  if sw == 2: # switch to print the cluster velocity
   print(t * 1E+15, v) # output time [fs] and electron cluster velocity [m / s]

  if sw == 3: # switch to print the cluster acceleration
   a = forceCompute(t, v, x, T) # cluster acceleration
   print(t * 1E+15, a) # output time [fs] and electron cluster acceleration [m / s2]

  if sw == 4: # switch to print the cluster temperature
   S = S + v ** 2 * h # calculating cluster temperature integral
   te = me * b * S / e # electron cluster temperature
   print(t * 1E+15, te) # output time [fs] and electron cluster temperature [ev]

  m1 = h * forceCompute(t, v, x, T) # velocity incriment in time t
  k1 = h * v # coordinate incriment in time t

  m2 = h * forceCompute(t + h / 2, v + m1 / 2, x + k1 / 2, T) # velocity incriment in time t + h/2
  k2 = h * (v + m1 / 2) # coordinate incriment in time t + h/2

  m3 = h * forceCompute(t + h / 2, v + m2 / 2, x + k2 / 2, T) # next velocity incriment in time t + h/2
  k3 = h * (v + m2 / 2) # next coordinate incriment in time t + h/2

  m4 = h * forceCompute(t + h, v + m3, x + k3, T) # next velocity incriment in time t
  k4 = h * (v + m3) # next coordinate incriment in time t

  t += h # next time
  v += (m1 + 2 * m2 + 2 * m3 + m4) / 6 # new velocity
  x += (k1 + 2 * k2 + 2 * k3 + k4) / 6 # new coordinate

 return

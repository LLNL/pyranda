# pyranda
A Python driven, Fortran powered Finite Difference solver for arbitrary hyperbolic PDE systems.  This is the mini-app for the Miranda code.


## Example use- solve the 1D advection equation
The one-dimensional advection equation is written as:

![Advection](http://mathurl.com/y7qnvzeg.png)

`from pyranda import pyrandaSim`

Define the domain/mesh

`domain = "xdom = (0.0 , 1.0 , 100 )"`

### Initialize a simulation object on a mesh
`pysim = pyrandaSim('advection',domain)`

### Define the equations of motion
`pysim.EOM(" ddt(:phi:)  =  -:c: * ddx(:phi:) ")`

#### Initialize variables
`pysim.setIC(":phi: = 1.0 + 0.1 * exp( -(abs(meshx-.5)/.1 )**2 )")`

### Integrate in time
`dt = .001`  
`time = 0.0`  
`while time < .1:`    
&nbsp;&nbsp;&nbsp;`time = pysim.rk4(time,dt)`  

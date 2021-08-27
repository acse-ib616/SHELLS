# 2D Finite Element Solver with Load Reconstruction

Foobar is a Python library for dealing with word pluralization.

The directory `MATLAB_elements` contains 4 scripts for the truss, CST, LST and TR elements programmed in MATLAB. It is not  the tool itself, but it has been useful for developing and checking implementations in a friendlier environment thant C++.

The directory `Script` contains the main script together with MATLAB functions and the `.cpp` files for the K and Jacobian sparse assemblers.

## Installation

Compile the `.cpp` files with the command:
```
mex .cpp
```

## Usage
Edit the `input.txt` file to vary the parametres of the model and to perform different operations:
```
Ly
nx
ny
E
t/r
nu
q
weight
uniform or random
cantilever or ss
CST/LST/Truss
stiffness
U
optimisation
tol
max it
lambda
plot
amp
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)

# Content
Python Implementation of the method presented in 
Im*proving a Constraint Programming Approach for Parameter  Estimation*.
Bertrand Neveu, Martin de la Gorce and Gilles Trombettoni,
27th IEEE International Conference on Tools with Artificial Intelligence 
Nov 9-11 2015, Vietri sul Mare, Italy 
Implemented By Martin de La Gorce , novembre 2015
Ecole des Ponts et Chaussees

#Installation
The code should run without the need to install dependencies. It uses numpy and matlplotlib.
You can run the line fitting example using from the command line 
using *python lineFitting.py*

# Example

We look for the paramters *a* and *b* such that 
|ax+b-y|<tolerance for at least three points (Q=3)
with tolerance=0.01

Points: 
![image0.png](https://bitbucket.org/repo/dEgXGz/images/274169658-image0.png)

Image of the search space after convergence:
![image62.png](https://bitbucket.org/repo/dEgXGz/images/3579631749-image62.png)
green rectangles either contain a solution or are too small to be refined
This is a tool program for aerospike rocket engine. 

Start the program with "Interface.py"

1. Step
The gas dynamic and thermal dynamic parameter can be inputted through RPA (Rocket propulsion analysis) or 
CEA (NASA Chemical Equilibrium with Applications) file. And for every new calculation a json will be created to save all inputted and calculated data.

2. Step
With the inputted gas dynamic, thermal dynamic parameter and neccesary geometrical parameter, 
the contour of aerospike thrust chamber can be generated, and the coordinates of the points on contour will be saved in the json file. 
or
You can input the existing contour data, which is built through the coordinates of points on the contour. 

3. Step
Input the wanted wall temperature, then the heat transfer around the aerospike thrust chamber can be calculated and be saved in the json file.

4. Step
You can check the contour of the generated aerospike thrust chamber, and you can get key parameters of the combustion chamber,
such as volume and characteristic length. 
And
You can check the heat transfer around the aerospike thrust chamber. And the total heat transfer will be shown in "kilowatt".

All the data can be outputted as txt file.

There are existing RPA and CEA file for testing.

All comtutational function in "Aerospike.py" can be also accessed for other uses.
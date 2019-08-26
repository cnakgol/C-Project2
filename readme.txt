This program finds the largest and the smallest eigenvalues (in norm) and the corresponding eigenvectors of a given real nxn matrix.
***To run the program, you should open the command line from your search button. It is better you change the direction to the program file.
Ex: This is written on your cmd:  C:\Users\Ceyda>
You should write "cd" and then write the file name containing the program: cd C:\Users\Ceyda\Desktop\CeydaNurAkgol
***After changing the location, you should call the program ".exe" and enter the inputs.
***3 inputs are required: 1-Input file name, 2-Tolerance, 3-Output file name
***Between each input, one space is required.
Ex: C:\Users\Ceyda\Desktop\CeydaNurAkgol>source.exe A.txt 1e-6 b.txt

*The given matrix should be nxn and real. This matrice should have been written in text document.
Example: 
2.7383  -0.5011  0.8817
-0.3039  2.3639  1.4258
0.1285  0.0665  3.8978 
It should also be in the same folder with the program. The name of the file should be given as the first input.
*A tolerance is needed as the second input. If you have no idea, you can enter 1e-6.
*The solution will be written to a file. Its name is needed as the third input. This file is written in the format below:
Eigenvalue#1: 4.00
0.38
0.80
1.00
Eigenvalue#2: 2.00
0.78
1,00
-0.09

*If the matrix is singular or nearly singular (any number which its absolute value is smaller than tolerance will be converted to 0),
a warning will appear in the console and no solution file will be generated.
*If the eigenvalues/eigenvectors have imaginary part, a warning will appear in the console and no solution file will be generated.
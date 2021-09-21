# matrix_solver
Gaussian Elimination with Partial Pivoting and Back Substitution

- Program gets two command line arguments. First one is name of the file that contains matrix A. Second one is name of the file that contains vector b.
- Program prints out the solution vector and saves it to a txt file called "x.txt".
- If matrix A is singular, program prints out that matrix is singular and quits.
- If matrix A is 2x2, program also calculates the condition numbers of 1 and infinity.

If condition number is too high, small changes in vector b will cause dramatic change in solutions x.

Example:

matrix A :  1.000 1.000        b1:   2.000     b2:  2.000
            1.000 1.001              2.000          2.001

Condition number of 1 of the matrix is 4004.
Condition number of infinity of the matrix is 4004.


Ax = b1  ---->  x1 = 2         Ax = b2  ---->   x2 = 1
                     0                               1

Conclusion:

0.001 change in vector b, causes solutions to change by 1.  (10^3 times change)

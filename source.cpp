#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#define Emach 1e-8  // arbitrarily defined a machine number 
using namespace std;

void exchangeRows(int firstRow, int secondRow, double **A, double b[],int dim); //row exchange function, used in partial pivoting
int upperTriangular(double **A, double b[],int dim);  // takes matrix and makes it upper triangular
void backSubstitution(double **A, double b[],double x[],int dim);  // using back substitution method, solves an upper triangular matrix
void solveMatrix(double **A, double b[],int dim);  // general function for solving matrix includes first 3 functions
void normCalc(double **A,double& absRowSum,double& absColumnSum); // calculates norm of the 2x2 matrices for condition numbers


int main (int argc,char *argv[]) {   // command line arguments given
  	
	ifstream file;   // file for input
	file.open(argv[1]); // opening file for taking matrix A input. First command line argument goes there
	string temp;
	int singularityCheck,dim = 0; // defining variables for checking singularity and taking dimension of matrix
	
	// reading matrix and vector from txt files
	
	// taking size of the matrix A
	while(getline(file,temp)){ 
		file >> temp;  // cursor is going through the file
		dim++;  // dim counter counts the number of loops so that matrix dimension is calculated
	}
	/////////////////////////////
	
	file.clear(); // cursor is at the end of the file
	file.seekg(0,ios::beg); // setting cursor to beginning of the file
	
	// creating dynamically allocated 2-d array matrix A
	double** matrixA = new double*[dim]; // firstly 1-d dynamic array is created
	
	for (int i=0; i<dim; i++) {
		matrixA[i] = new double[dim];  // second dimension of the dynamic array is created
	}
	/////////////////////////////
	
	// copying entries from txt file to the array
	for(int i=0; i<dim; i++){
		for(int j=0; j<dim; j++){
			file >> matrixA[i][j];  //looping through file and copying matrix entries to our dynamic array
		}
	}
	file.close(); // closing file
	/////////////////////////////    
	
	// creating dynamically allocated vector b
	double* b = new double[dim];   // dynamically allocating for array b
	
	file.open(argv[2]);  // opening the second file. Second command line argument goes here.
	for(int i=0; i<dim; i++){
		file >> b[i];    // copying contents of the file to array b
	}
	file.close(); // closing file
	///////////////////////////////////////////////////////////
	
	// solution process
	if (dim == 2){  // if matrix is 2x2, we will also calculate condition numbers 
		
		double determinant = (matrixA[0][0] * matrixA[1][1]) - (matrixA[1][0] * matrixA[0][1]); // calculating determinant of the 2x2 matrix
		
		if (determinant < Emach && determinant > -1*Emach){  // if determinant is smaller than our Emach (nearly 0)
			cout << "Condition number of this singular 2x2 matrix is infinity by definition.";  // print out that matrix is singular and by convention its condition number is infinity
			return 1;
		}
		
		double** inverseA = new double*[2]; // if matrix is nonsingular we declare dynamically allocated 2-d matrix for A inverse
		
		for (int i=0; i<2; i++) {
			inverseA[i] = new double[2]; // creating second dimension of the matrix
		}
		
		inverseA[0][0] = matrixA[1][1]/determinant;    // inverse of 2x2 matrices is trivial. Just change of 1st and 4th entries and sign change of other entries. After that divide by determinant
		inverseA[0][1] = -1*matrixA[0][1]/determinant;
		inverseA[1][0] = -1*matrixA[1][0]/determinant;
		inverseA[1][1] = matrixA[0][0]/determinant;
		
		
		double temp1,temp2,rowSumA,columnSumA,rowSumInverse,columnSumInverse; // declaring necessary variables to store information for calculating condition numbers
		
		normCalc(matrixA,temp1,temp2);  // calculating 1 and infinity norm of matrix A
		rowSumA = temp1;  // temp1 is assigned to rowSumA
		columnSumA = temp2; // temp2 is assigned to columnSumA
		
		normCalc(inverseA,temp1,temp2); // calculating 1 and infinity norm of inverse matrix A
		rowSumInverse = temp1; // temp1 is assigned to rowSumInverse
		columnSumInverse = temp2; // temp2 is assigned to columnSumInverse
		
		cout << "Condition number of 1 of the matrix is " << columnSumA * columnSumInverse << endl;  // printing out condition numbers
		cout << "Condition number of infinity of the matrix is " << rowSumA * rowSumInverse << endl;
		
		
		for (int i= 0; i < 2; i++) { // deleting inverse matrix a from memory
			delete[] inverseA[i]; // deleting second dimension entries
		}
		delete[] inverseA; // deleting matrix
		
		
	}
	
	singularityCheck = upperTriangular(matrixA,b,dim); // making matrix upper triangular and return value is assigned to singularityCheck
	
	if(singularityCheck){ // if function returns 1 and singularity Check becomes 1
		cout << "This matrix is singular!!!";  // printing out matrix is singular
	}
	else{
		solveMatrix(matrixA,b,dim);  // solution process of matrix A
	}
	///////////////////////////////////////////////////////////
	
	// deleting 2-d array matrix A
	for (int i= 0; i < dim; i++) { 
		delete[] matrixA[i]; // deleting second dimension of matrix A
	}
	delete[] matrixA; // deleting matrix A
	
	// deleting vector b
	delete[] b;
	
	return 0;
}

void exchangeRows(int firstRow, int secondRow, double **A, double b[],int dim){
	
	double temp;
	for(int i=0; i<dim; i++){
		// exchanging entries of matrix A
		temp = A[firstRow][i];
		A[firstRow][i]	= A[secondRow][i];
		A[secondRow][i] = temp;
	}
	// exchanging entries of vector b
	temp = b[firstRow];
	b[firstRow] = b[secondRow];
	b[secondRow] = temp;
}

int upperTriangular(double **A, double b[],int dim){
	
	double temp;
	
	// upper trianglization process
	for(int k=0; k<dim-1; k++){ // implementing row operations to all dimensions
		
		double largestPivot = k;    // we need to find largest absolute value of the column and select it as pivot. (partial pivoting)
		for(int i=k; i<dim; i++){   // looping through column
			if (abs(A[i][k]) > abs(A[k][k])){  // if we find a larger entry 
				largestPivot = i;  // we set largest pivot to that entry
			}
		}
		exchangeRows(k,largestPivot,A,b,dim);  // exchanging rows according to the partial pivoting
		
		for(int j=k+1; j<dim; j++){ // implement row operation to all rows
			
			temp = ((-1*A[j][k])/(A[k][k]));  // temp is set to the value that multiplies the row and adds it to another
			
			for(int i=0; i<dim; i++){  // implement row operation to all column entries
				A[j][i] = A[k][i] * temp + A[j][i]; // implementing row operation to all columns of matrix A
			}
			
			b[j] = ((b[k]) * temp) + b[j]; // implementing same row operation to vector b
			
		}
		
	}
	
	// singularity check
	for(int i=0; i<dim; i++){ // looping through all pivot points
		if(A[i][i] < Emach && A[i][i] > -1 * Emach)  // if any pivots is smaller than Emach meaning nearly 0
			return 1; // we return 1 to indicate that matrix is singular
	}
	return 0; // otherwise return 0
	
}

void backSubstitution(double **A, double b[],double x[],int dim){
	
	// after calculating upper triangular, we implement back substitution to find solutions
	double knownSolutions = 0; // this variable holds the previously find solutions x and their coefficients for calculating the current x
	
	for(int i=dim-1; i>=0; i--){ // looping backwards this time
		for(int j=dim-1; j>i; j--){
			knownSolutions += A[i][j] * x[j]; // multiplying known solutions with their coefficients and adding them up
		}
		x[i] = ((b[i]) - knownSolutions) / (A[i][i]);  // to find current x value we have substituted known solutions and leave x alone in the equation
		knownSolutions = 0; // after the solution step, we need to set this variable to 0 in order to calculate the next x value
	}
	
}

void solveMatrix(double **A, double b[],int dim) {
	
	double x[dim]; // creating array for storing solutions
	
	
	// back substitution process
	backSubstitution(A,b,x,dim); // implement back substitution (upper triangular process already made)
	
	
	// printing x vector
	for(int i=0; i<dim; i++){  
		cout << x[i] << endl; // printing out all solutions in correct order
	}
	
	// writing the solutions to x.txt
	ofstream myfile; 
  	myfile.open ("x.txt");  // creating a file named "x.txt" to write the solutions
  	for(int i=0;i<dim;i++){
		myfile << x[i] << endl; // writing the solutions
	} 
	myfile.close();  // closing file
}

void normCalc(double **A,double& absRowSum,double& absColumnSum){
	
	// to calculate 1 and inf norms of matrices, we need to find absolute row sum and absolute column sum
	// we need to calculate 2x2 matrix norms only. So we have 2 rows and 2 columns to compare
	
	if (abs(A[0][0])+abs(A[0][1]) < abs(A[1][0])+abs(A[1][1])){ 
		absRowSum = abs(A[1][0])+abs(A[1][1]); // row sum is set according to the larger row
	}
	else{
		absRowSum = abs(A[0][0])+abs(A[0][1]); // row sum is set according to the larger row
	}
	
	if (abs(A[0][0])+abs(A[1][0]) < abs(A[0][1])+abs(A[1][1])){
		absColumnSum = abs(A[0][1])+abs(A[1][1]); // column sum is set according to the larger row
	}
	else{
		absColumnSum = abs(A[0][0])+abs(A[1][0]); // column sum is set according to the larger row
	}
}
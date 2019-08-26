#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
using namespace std;

class Matrix{
	float **mat;// creates a variable to hold the content of the matrix
	public:
		int size;// creates a variable to hold the size of the matrix
		Matrix(){size=0;}; // creates default constructor
		Matrix(const Matrix&);// creates copy constructor
		void findSize(string);// reads the matrix from a file and finds the number of lines i.e. size for a square matrix
		void read(string);// reads the matrix from a file and writes it to a matrix
		void multiplication(float*, float*);// multiplies a matrix with a vector
		void partialPivoting(float*,float);// does partial pivoting
		int backwardSubstitution(float*,float*,float);// makes backward substitution
		int isDiagonal();// checks whether the matrix is diagonal or not
		float kthMember(int, int);// finds the element of the matrix based on the given coordinates
		Matrix& operator=(const Matrix&);// equals matrices
		~Matrix();// creates destructor
};

Matrix::Matrix(const Matrix &param){
	size=param.size;// equals size
	mat=new float* [size];
	for(int i=0; i<size; i++){
		mat[i]=new float[size];
	}// creates dynamic allocation for the new matrix
	for(int i=0; i<size; i++){
		for(int j=0; j<size; j++){
			mat[i][j]=param.mat[i][j];
		}
	}// equals matrices
}

void Matrix::findSize(string inputfile){
	string line;
	ifstream matrixfile(inputfile.c_str());// creates objects of ifstream class to make operations with the input file, i.e. reading from the input file

	if(matrixfile.is_open()){
		while(getline(matrixfile,line)){
			size++;
		}// counts how many lines in the file, which indicates the size of matrix
		matrixfile.close();// closes the file
	} 
	else cout << "Unable to open the input file." << "\n"; // warns if the input file couldn't be opened.
	//cout << n; to control whether it counts the lines correctly or not
}

void Matrix::read(string inputfile){
	mat=new float* [size];
	for(int i=0; i<size; i++){
		mat[i]=new float[size];
	} // creates dynamic allocation for the matrix A
	string line;
	ifstream matrixfile (inputfile.c_str());// creates objects of ifstream class to make operations with the input file, i.e. reading from the input file
	
	if(matrixfile.is_open()){
		for(int i=0; i<size; i++){
			for(int j=0; j<size; j++){
				matrixfile >> mat[i][j];
			}
		}
	}// reads the input and write it to the 2d array
	else cout << "Unable to open the input file." << "\n"; // warns if the input file couldn't be opened.
}

void Matrix::multiplication(float* vec1, float* vec2){
	for(int i=0; i<size; i++){
		for(int j=0; j<size; j++){
			vec2[i]+=mat[i][j]*vec1[j];
		}
	 }
}// multiplies a matrix nxn with a vector nx1 and returns the solution vector nx1

void Matrix::partialPivoting(float* vec,float precision){
	float max;// creates a variable to keep the maximum element of the every column for partial pivoting
	int maxindice;// creates a variable to keep the indice of the maximum element of the every column for partial pivoting
	
	for(int j=0; j<size-1; j++){// loop in columns except the last one which doesn't require pivoting
		max=abs(mat[j][j]);// set the maximum to the diagonal of the current column
		maxindice=j;// set the maximum indice to the current row/column
		
		for(int i=j+1; i<size; i++){// loop in rows except the first one which is already kept with the max
			if(abs(mat[i][j])>max){// compares whether the current entry is bigger than the maximum value of the column
				max=abs(mat[i][j]);// sets maximum to the current entry if it is bigger than the current maximum value
				maxindice=i;// keeps the indice/row number of the new maximum value
			}
		}
		if(maxindice!=j){// controls whether the maximum is changed
			swap(mat[j],mat[maxindice]);// if changed, changes the rows of the matrix A
			swap(vec[maxindice],vec[j]);// if changed, changes the rows of the matrix b
			//maxindice=j+1;// set the maximum indice for the next column
		}
		
		float mult;// creates a variable to keep the scaling factor/multiplier for Gaussian elimination
		
		for(int i=j+1; i<size; i++){
			if(abs(mat[j][j])<precision){
				break;// exits the loop whether the matrix is singular i.e. after the partial pivoting, any diagonal entry is zero
			}
			else mult=mat[i][j]/mat[j][j];// if not, calculates the scaling factor/multiplier for Gaussian elimination
				
			for(int k=0; k<size; k++){
				mat[i][k]=mat[i][k]-mult*mat[j][k];// set A's entries to zero
			}
				
			vec[i]=vec[i]-mult*vec[j];// updates also b's entries
		}	
	}
}

int Matrix::backwardSubstitution(float* vec1, float* vec2, float precision){
	for(int j=size-1; j>=0; j--){
		if(abs(mat[j][j])<precision){
			cout << "The matrix is singular.";
			return 0;// warns if the matrix is singular and quits
		}
		vec2[j]=vec1[j]/mat[j][j];
		
		for(int i=0; i<j; i++){
			vec1[i]=vec1[i]-mat[i][j]*vec2[j];
		}// solves the upper triangular matrix and writes the solution to a given vector 
	}
	return 1;
}

int Matrix::isDiagonal(){
	for(int i=0; i<size; i++){
		for(int j=0; j<size; j++){
			if(i!=j){
				if(fabs(mat[i][j])>=1e-6) return 0;
			}// checks whether the non-diagonal elements are non-zero (or close to non-zero)
		}
	}
	return 1;
}

Matrix& Matrix::operator=(const Matrix& param){
	if(this!=&param){
		size=param.size;// equals size
		mat=new float* [size];
		for(int i=0; i<size; i++){
			mat[i]=new float[size];
		} // creates dynamic allocation for the new matrix
		for(int i=0; i<size; i++){
			for(int j=0; j<size; j++){
				mat[i][j]=param.mat[i][j];
			}
		}// equals matrices
	}
	return *this;
}

float Matrix::kthMember(int a, int b){
	return mat[a][b];// finds the element of the matrix based on the given coordinates
}

Matrix::~Matrix(){
	for(int i=0; i<size; i++){
			delete[] mat[i];
		}
	delete[] mat;// deletes the matrix
}

class Vector{
	public:
		int size;// creates a variable to hold the size of the vector
		float *vec;// creates a variable to hold the content of the vector
		float norm;// creates a variable to hold the infinity norm of the vector
		void getSize(const Matrix&);// gets the size of a matrix to construct a same sized vector
		void generateUnitVector();// generates unit vector with the first element 1 and all the others 0
		void generateNullVector();// generates null vector with all the elements 0
		void findNorm();// finds the infinity norm of a vector, i.e. the absolute maximum element
		void normalization();// normalizes the vector by dividing all the elements to the norm
		Vector& operator=(const Vector&);// overloads the assignment operator
		int isEqual(const Vector&, float);// checks whether the two vectors are the same or not
		int isAbsoluteInverse(const Vector&, float);// checks whether the absolute values of the elements of two vectors are the same or not
		/*this function is written in order to avoid problems due to negative matrices' non-convergence problems 
		i.e. [-1 -1]T and [1 1]T belong to the same group of eigenvectors but cannot converge to each other*/
		~Vector(){delete vec;}// deletes the vector
};

void Vector::getSize(const Matrix &param){
	size=param.size;// assigns the size of the matrix to the vector
}

void Vector::generateUnitVector(){
	vec=new float[size];// creates a dynamically allocated vector
	vec[0]=1;// equals the first element to 1
	for(int i=1; i<size; i++){
		vec[i]=0;
	}// generates unit vector with the elements 0
}

void Vector::generateNullVector(){
	vec=new float[size];// creates a dynamically allocated vector
	for(int i=0; i<size; i++){
		vec[i]=0;// equals all the elements to 0
	}
}

void Vector::findNorm(){
	float max=fabs(vec[0]);// creates a variable to hold the element of maximum amplitude and equals it to the first element to compare with others
	for(int i=1; i<size; i++){
		max=fmax(fabs(vec[i]),max);// checks whether any other element is greater than the max value and if so, assign it to max value
	}
	norm=max;// returns max value as infinity norm
}

void Vector::normalization(){
	for(int i=0; i<size; i++){
		vec[i]=vec[i]/norm;
	}// divides all elements to norm in order to normalize the vector
}

Vector& Vector::operator=(const Vector& param){
	if(this!=&param){
		size=param.size;// equals sizes
		delete vec;
		vec=new float[size];
		for(int i=0; i<size; i++){
			vec[i]=param.vec[i];
		}// equals vectors
		norm=param.norm;// equals norms
	}
	return *this;
}

int Vector::isEqual(const Vector& param, float precision){
	for(int i=0; i<size; i++){
		if(fabs(param.vec[i]-vec[i])>precision){
			return 0;
		}
	}
	return 1;
}// checks whether the two vectors are the same or not

int Vector::isAbsoluteInverse(const Vector& param, float precision){
	for(int i=0; i<size; i++){
		if(fabs(vec[i]+param.vec[i])>=precision){
			return 1;
		}
	}
	return 0;
}// checks whether the absolute values of the elements of two vectors are the same or not


void solutionFile(string outputfile, Vector& vec1, Vector& vec2, float val1, float val2){
	ofstream solutionfile(outputfile.c_str());// creates objects of ofstream class to make operations with output file, i.e. writing to it

	if(solutionfile.is_open()){
		solutionfile << "Eigenvalue#1: " << val1 << endl;
		for(int i=0; i<vec1.size; i++){
			solutionfile << vec1.vec[i] << endl;
		}// writes the eigenvalue of the largest amplitude and the corresponding eigenvector to the solution file
		solutionfile << "Eigenvalue#2: " << val2 << endl;
		for(int i=0; i<vec2.size; i++){
			solutionfile << vec2.vec[i] << endl;
		}// writes the eigenvalue of the smallest amplitude and the corresponding eigenvector to the solution file
		solutionfile.close();// closes the solution file
	}
	else cout << "Unable to open the solution file." << "\n"; // warns if the solution file couldn't be opened.
}

int main(int argc, char *argv[]){
	string inputfile="A.txt";
	float precision=1e-6;
	string outputfile="x.txt";
	
	if(argc>=1){
		inputfile=argv[1];// receives the first entry of the user's input as the name of the input file
	}
	if(argc>=2){
		precision=atof(argv[2]);// receives the second entry of the user's input as the precision
	}
	if(argc>=3){
		outputfile=argv[3];// receives the third entry of the user's input as the name of the output file
	}
	
	Matrix A;// creates the main matrix
	A.findSize(inputfile);// finds its size
	A.read(inputfile);// reads the content
	
	Vector b1;// creates the vector x(k-1) for the normalized power iteration
	b1.getSize(A);// gets its size from the main matrix
	b1.generateUnitVector();// generates a non-zero vector i.e. unit vector
	
	Vector b2;// creates the vector x(k-1) for the inverse power iteration
	b2.getSize(A);// gets its size from the main matrix
	b2.generateUnitVector();// generates a non-zero vector i.e. unit vector
	
	Vector x1;// creates a vector for the eigenvector corresponding to the largest eigenvalue
	x1.getSize(A);// gets its size from the main matrix
	x1.generateNullVector();// generates a zero vector
	
	Vector x2;// creates a vector for the eigenvector corresponding to the smallest eigenvalue
	x2.getSize(A);// gets its size from the main matrix
	x2.generateNullVector();// generates a zero vector
	
	float eigenval1=0;// creates a variable for the largest eigenvalue
	float eigenval2=0;// creates a variable for the smalleest eigenvalue
	
	if(A.isDiagonal()){// calculates all the values if the matrix is diagonal
		float max=fabs(A.kthMember(0,0));// creates a variable for the maximum value in the matrix and equals it to the first element
		float min=fabs(A.kthMember(0,0));// creates a variable for the minimum value in the matrix and equals it to the first element
		int indice1=0;// creates a variable to keep the max value's place and equals it to 0
		int indice2=0;// creates a variable to keep the min value's place and equals it to 0
		for(int i=0; i<A.size; i++){
			if(fabs(A.kthMember(i,i))>max){// checks if the current value is greater than the max value
				max=fabs(A.kthMember(i,i));// assigns it to max value
				indice1=i;// assign its indice to the max indice
			}
			if(fabs(A.kthMember(i,i))<min){// checks if the current value is smaller than the min value
				min=fabs(A.kthMember(i,i));// assigns it to min value
				indice2=i;// assign its indice to the min indice
			}
		}		
		x1.vec[indice1]=1;// turns x1 into the corresponding unit vector
		x2.vec[indice2]=1;// turns x2 into the corresponding unit vector
		eigenval1= max;// assigns max value to the largest eigenvalue
		eigenval2= min;// assigns min value to the smallest eigenvalue
		solutionFile(outputfile, x1, x2, eigenval1, eigenval2);// writes the solutions to the output file
		
		return 0;// ends the program
	}

	int maxIteration=0;// creates a variable to keep iteration number

	while(1){// the normalized power iteration
		maxIteration++;// counts the iterations
		if(maxIteration>=1000){
			cout << "Error: The values don't converge.";
			return 0;
		}// ends the program if the iteration number is greater than 500 i.e. values don't seem to be converged
		
		A.multiplication(b1.vec,x1.vec);// multiplies A by b1, writes the result to x1
		x1.findNorm();// finds norm of x1
		eigenval1=x1.norm;// assigns the norm to the largest eigenvalue
		x1.normalization();// normalizes x1
		
		if(b1.isEqual(x1,precision)==1){
			break;
		}// ends the loop if the eigenvector converges
		if(b1.isAbsoluteInverse(x1,precision)==0){
			eigenval1=-eigenval1;
			break;
		}// ends the loop if the eigenvector is changing each time to its additive inverse (explained in the isAbsoluteInverse function above)
		b1.operator=(x1);// assigns x1 to b1 to continue new iterations
		x1.generateNullVector();// empties x1
	}
	
	Matrix M(A);// creates a new matrix to keep A (because it changes after partial pivoting but the unchanged version is required for the inverse power iteration)
	
	Vector c;// creates a new vector to keep b2 (because it changes after partial pivoting but the unchanged version is required for the inverse power iteration)
	c.getSize(A);// gets its size from the main matrix
	c.generateUnitVector();// generates a zero vector
	
	maxIteration=0;// equals iteration number to zero in order to count for the inverse power iteration

	while(1){// the inverse normalized power iteration
		maxIteration++;// counts the iterations
		if(maxIteration>=1000){
			cout << "Error: The values don't converge.";
			return 0;
		}// ends the program if the iteration number is greater than 500 i.e. values don't seem to be converged
		
		A.partialPivoting(b2.vec, precision);
		if(A.backwardSubstitution(b2.vec,x2.vec,precision)==0){
	 		return 0;
		}// solves Ax2=b2 after partial pivoting and backward subsitution
		
	 	x2.findNorm();// finds norm of x2
	 	eigenval2=(1/x2.norm);// assigns the norm's multiplication inverse to the smallest eigenvalue
	 	x2.normalization();// normalizes x2
	 	
	 	if(c.isEqual(x2,precision)==1){
			break;
		}// ends the loop if the eigenvector converges
		if(c.isAbsoluteInverse(x2,precision)==0){
			eigenval2=-eigenval2;
			break;
		}// ends the loop if the eigenvector converges
		b2.operator=(x2);// assigns x2 to b2 to continue new iterations
		c.operator=(b2);// assigns b2 to c to continue new iterations
		A.operator=(M);// updates A from the unchanged M matrix
	}
		
	solutionFile(outputfile, x1, x2, eigenval1, eigenval2);// writes solutions to a file which its name is given
}

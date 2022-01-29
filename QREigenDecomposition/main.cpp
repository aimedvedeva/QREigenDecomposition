#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>

using namespace std;

ifstream fileOut;
ofstream file;

static double *matrix;
static int numOfVariables;
static int numOfEquations;

static int numOfGivensMultIter = 0;
static int numOfQRSteps = 0;

void GetMatrix(int numOfVar, int numOfEq, double *matrix/*out*/) {
	for (int i = 0; i < numOfEq; i++) {
		for (int j = 0; j < numOfVar; j++) {
			fileOut >> matrix[i * numOfVar + j];
		}
	}
}

void GetMatrixSize(int *numOfVar/*out*/, int *numOfEq/*out*/) {
	fileOut >> *numOfEq >> *numOfVar;
	return;
}

void PrintMatrix(double *matrix, int numOfVariables, int numOfEquations) {
	file.precision(6);

	for (int i = 0; i < numOfEquations; i++) {
		for (int j = 0; j < numOfVariables; j++) {
			file << matrix[i * numOfVariables + j] << "     ";
		}
		file << endl;
	}
	file << "--------------------------------------";
	file << endl;
}

double *MultMatrixAB(double *matrixA, int numOfVariablesA, int numOfEquationsA, double *matrixB, int numOfVariablesB, int numOfEquationsB) {
	if (numOfVariablesA != numOfEquationsB) {
		return NULL;
	}

	double *newMatrix = new double[numOfEquationsA * numOfVariablesB];

	for (int i = 0; i < numOfEquationsA; i++) {
		for (int j = 0; j < numOfVariablesB; j++) {
			double sum = 0;
			for (int k = 0; k < numOfVariablesA; k++) {
				sum += matrixA[i * numOfVariablesA + k] * matrixB[k * numOfVariablesB + j];
			}
			newMatrix[i * numOfVariablesB + j] = sum;
		}
	}
	
	return newMatrix;
}

double EuclideanNorm(double *vector, int lenght) {
	double sum = 0;
	for (int i = 0; i < lenght; i++) {
		sum += powl(vector[i], 2);
    }
	sum = sqrt(sum);
	return sum;
}

double *TransposeMatrix(double *matrix, int numOfColumns, int numOfRows)
{
	double *newMatrix = new double[numOfColumns * numOfRows];

	for (int i = 0; i < numOfRows; i++) {
		for (int j = 0; j < numOfColumns; j++) {
			newMatrix[j * numOfColumns + i] = matrix[i * numOfColumns + j];
		}
	}

	return newMatrix;
}

void SubtractMatrix(double *matrixA, double *matrixB, int numOfColumns, int numOfRows) {
	
	for (int i = 0; i < numOfRows; i++) {
		for (int j = 0; j < numOfColumns; j++) {
			matrixA[i * numOfColumns + j] -= matrixB[i * numOfColumns + j];
		}
	}
}

void CopyMatrix(double *to, double *from, int numOfRows, int numOfColumns, int fromColumn, int fromRow) {
	for (int i = fromRow; i < numOfRows; i++) {
		for (int j = fromColumn; j < numOfColumns; j++) {
			to[i * numOfColumns + j] = from[i * numOfColumns + j];
		}
	}
}

double *NormalForHausholderDec(double *matrixA, int numOfColumns, int numOfEquations, int step) {
	double *normal = new double[numOfEquations - step - 1];

	//copy the column we need
	for (int i = step + 1; i < numOfEquations; i++) {
		normal[i - step - 1] = matrixA[i * numOfColumns + step];
	}

	double factor = EuclideanNorm(normal, numOfEquations - step - 1);

	if (normal[0] >= 0) {
		factor = -factor;
	}

	normal[0] -= factor;//simulating the subtraction of the first vector of the natural basis * factor

	return normal;
}

void InitByZero(double *data, int size) {
	for (int i = 0; i < size; i++) {
		data[i] = 0;
	}
}

double *GetIdentityMatrix(int numOfColumns, int numOfRows) {
	double *identityMatrix = new double[numOfColumns * numOfRows];
    
	
	InitByZero(identityMatrix, numOfColumns * numOfRows);
	for (int i = 0; i < numOfRows; i++) {
		for (int j = 0; j < numOfColumns; j++) {
			if (i == j) {
				identityMatrix[i * numOfColumns + j] = 1;
			}
		}
	}
	return identityMatrix;
}

void MultMatrixByNum(double *matrix, int numOfColumns, int numOfRows, double num) {
	for (int i = 0; i < numOfRows; i++) {
		for (int j = 0; j < numOfColumns; j++) {
			matrix[i * numOfColumns + j] *= num;
		}
	}
}

double *ReflectionMatrix(double *normal, int lenght) {
	double norm = EuclideanNorm(normal, lenght);

	double *reflectionMatrix = GetIdentityMatrix(lenght, lenght);

	double *transposeNormal = TransposeMatrix(normal, 1, lenght);

	double *mult = MultMatrixAB(normal, 1, lenght, transposeNormal, lenght, 1);//w * wT
	delete[] transposeNormal;

	MultMatrixByNum(mult, lenght, lenght, (2 / powl(norm, 2)));//2/norm^2 * (w * wT)

	SubtractMatrix(reflectionMatrix, mult, lenght, lenght);

	free(mult);

	return reflectionMatrix;
}

double *GetMatrixForEachStep(double *reflectionMatrix, int numOfColumns, int numOfRows, int step) {
	double *resultMatrix = GetIdentityMatrix(numOfColumns, numOfRows);

	for (int i = step; i < numOfRows; i++) {
		for (int j = step; j < numOfColumns; j++) {
			resultMatrix[i * numOfColumns + j] = reflectionMatrix[(i - step) * (numOfColumns - step) + (j - step)];
		}
	}

	return resultMatrix;
}

void HouseholderTransformation() {

	//loop through columns
	for (int i = 0; i < numOfVariables - 2; i++) {
		double *normal = NormalForHausholderDec(matrix, numOfVariables, numOfEquations, i);
		double *reflectionMatrix = ReflectionMatrix(normal, numOfEquations - i - 1);
		double *factorMatrix = GetMatrixForEachStep(reflectionMatrix, numOfVariables, numOfEquations, i + 1);

		double *reflectionMultMatrix = MultMatrixAB(factorMatrix, numOfVariables, numOfEquations, matrix, numOfVariables, numOfEquations); /* P * A = B */
		double *matrixMultReflection = MultMatrixAB(reflectionMultMatrix, numOfVariables, numOfEquations, factorMatrix, numOfVariables, numOfEquations); /* B * P */

        delete[] matrix;
		delete[] reflectionMultMatrix;

		matrix = matrixMultReflection;
		
		delete[] normal;
		delete[] reflectionMatrix;
	}
}


//-----------------------------------------
void CalculateRotationMatrixParameters(double *matrix, double *c /*cos(a)*/, double *s /*sin(a)*/, int j) {
	double t = matrix[j * numOfVariables + j] / matrix[(j + 1) * numOfVariables + j];
	*c = 1 / sqrt(1 + powl(t, 2));
	*s = t * (*c);
}

double *GivensTransformation() {
	double *G = GetIdentityMatrix(numOfVariables, numOfEquations); // G1^(-1) * G2^(-1) * ...

	//loop through columns
	for (int i = 0; i < numOfVariables - 1; i++) {
		double c, s;
		CalculateRotationMatrixParameters(matrix, &c, &s, i);
		double *rotationMatrix = GetIdentityMatrix(numOfVariables, numOfEquations);

		rotationMatrix[i * numOfVariables + i] = s;
		rotationMatrix[(i + 1) * numOfVariables + (i + 1)] = s;

		rotationMatrix[i * numOfVariables + (i + 1)] = c;
		rotationMatrix[(i + 1) * numOfVariables + i] = -c;

		double *newMatrix = MultMatrixAB(rotationMatrix, numOfVariables, numOfEquations, matrix, numOfVariables, numOfEquations);

		delete[] matrix;
		matrix = newMatrix;


		double *transposeRotMatrix = TransposeMatrix(rotationMatrix, numOfVariables, numOfEquations);
		double *tmpMatrix = MultMatrixAB(G, numOfVariables, numOfEquations, transposeRotMatrix, numOfVariables, numOfEquations);

		delete[] transposeRotMatrix;
		delete[] G;

		G = tmpMatrix;
		numOfGivensMultIter++;
    }
	return G;
}

bool IsTrueCondToStopIter(double error) {
	double *vector = new double[numOfVariables - 1];
	for (int i = 1; i < numOfVariables; i++) {
		vector[i - 1] = matrix[i * numOfVariables + (i - 1)];
	}
	double norm = EuclideanNorm(vector, numOfVariables - 1);

	delete[] vector;
	if (norm < error) {
		return true;
	}
	return false;
}

void QREigenDecomposition(double error) {
	bool isFirstIter = true;

	while (!(IsTrueCondToStopIter(error)) || isFirstIter) {
		double *G = GivensTransformation(); // получили R и G, A(k) = G * R
		
		double *tmpMatrix = MultMatrixAB(matrix, numOfVariables, numOfEquations, G, numOfVariables, numOfEquations); // A(k + 1) = R * G
		delete[] matrix;
		matrix = tmpMatrix;

		numOfQRSteps++;
		isFirstIter = false;
    }
}


//-----------------------------------------
double *NormalForReflectionMethod(double *matrixA, int numOfColumns, int numOfEquations, int step) {
	double *normal = new double[numOfEquations - step];

	// copy the column we need
	for (int i = step; i < numOfEquations; i++) {
		normal[i - step] = matrixA[i * numOfColumns + step];
	}

	double factor = EuclideanNorm(normal, numOfEquations - step);

	if (normal[0] >= 0) {
		factor = -factor;
	}

	normal[0] -= factor;//imitate the subtraction of the first vector of the natural basis * factor

	return normal;
}

double *BackElimination(double *matrix, double *constTermsVec, int numOfVariables, int numOfEquations) {
	double *solutions = new double[numOfVariables];
	int numOfSolutions = numOfVariables - 1;
	int curRow = numOfEquations - 1;
	int curColumn = numOfVariables - 1;

	while (curRow >= 0) {
		double sum = 0;
		for (int i = curColumn + 1; i < numOfVariables; i++) {
			sum += matrix[curRow * numOfVariables + i] * solutions[i];
		}
		solutions[numOfSolutions] = (constTermsVec[curRow] - sum) / matrix[curRow * numOfVariables + curColumn];
		curRow--;
		curColumn--;
		numOfSolutions--;
	}
	return solutions;
}

double *ReflectionMethod(double *matrixA, double *matrixB, int numOfVariables, int numOfEquations) {

	//loop through columns
	for (int i = 0; i < numOfVariables - 1; i++) {
		double *normal = NormalForReflectionMethod(matrixA, numOfVariables, numOfEquations, i);
		double *reflectionMatrix = ReflectionMatrix(normal, numOfEquations - i);
		double *factorMatrix = GetMatrixForEachStep(reflectionMatrix, numOfVariables, numOfEquations, i);
		double *newMatrixA = MultMatrixAB(factorMatrix, numOfVariables, numOfEquations, matrixA, numOfVariables, numOfEquations);
        double *newMatrixB = MultMatrixAB(factorMatrix, numOfVariables, numOfEquations, matrixB, 1, numOfEquations);

		delete[] matrixA;
		delete[] matrixB;

		matrixA = newMatrixA;
		matrixB = newMatrixB;

		delete[] normal;
		delete[] reflectionMatrix;
		delete[] factorMatrix;
	}
	double *X = BackElimination(matrixA, matrixB, numOfVariables, numOfEquations);
	delete[] matrixA;
	delete[] matrixB;
	return X;
}


//-----------------------------------------
double *CreateVectorWithSingleNorm(int lenght) {
	double *vector = new double[lenght];

	for (int i = 0; i < lenght; i++) {
		vector[i] = 1 / sqrt(lenght);
	}

	return vector;
}

double DotProduct(double *matrix, double *vector) {
	//(Ax,x)
	double *AX = MultMatrixAB(matrix, numOfVariables, numOfEquations, vector, 1, numOfEquations);

	double result = 0;

	for (int i = 0; i < numOfEquations; i++) {
		result += AX[i] * vector[i];
	}
	
	delete[] AX;

	return result;
}

double *RayleighQuotientIteration(double eigenNumber, double *matrix, double C) {
	double *X  = CreateVectorWithSingleNorm(numOfEquations);

	double normY = 0;

	while (normY < C) {
		double *AMinusEigenNumE = new double[numOfEquations * numOfVariables];
		CopyMatrix(AMinusEigenNumE, matrix, numOfEquations, numOfVariables, 0, 0);

		double *EigenNumE = GetIdentityMatrix(numOfVariables, numOfEquations);
		MultMatrixByNum(EigenNumE, numOfVariables, numOfEquations, eigenNumber);

		SubtractMatrix(AMinusEigenNumE, EigenNumE, numOfVariables, numOfEquations);
		delete[] EigenNumE;

		double *Y = ReflectionMethod(AMinusEigenNumE, X, numOfVariables, numOfEquations);

		normY = EuclideanNorm(Y, numOfEquations);

		X = Y;

        MultMatrixByNum(X, 1, numOfEquations, 1 / normY);

		eigenNumber = DotProduct(matrix, X);
    }
    
	return X;
}


int main(void) {
	fileOut.open("matr1.txt");

	GetMatrixSize(&numOfVariables, &numOfEquations);

	//calculate the matrix
	matrix = new double[numOfVariables * numOfEquations];
	GetMatrix(numOfVariables, numOfEquations, matrix);

	fileOut.close();

	//--------------------------------------------------

	//create a copy of the original matrix for reverse iterations
	double *copyMatrix = new double[numOfEquations * numOfVariables];
	CopyMatrix(copyMatrix, matrix, numOfEquations, numOfVariables, 0, 0);


	//reduce the original matrix to the Hessenberg form
	HouseholderTransformation();


	//find own numbers QR - method
	double error = 0.0001;
	QREigenDecomposition(error);

	file.open("result.txt");


	//print the resulting triangular matrix with similar eigenvalues
	PrintMatrix(matrix, numOfVariables, numOfEquations);

	double C = 50000;

	//in the loop by eigenvalues ​​we find the relative vectors -> print to a file
	for (int i = 0; i < numOfEquations; i++) {
		double *eigenVector = RayleighQuotientIteration(matrix[i * numOfVariables + i], copyMatrix, C);
		PrintMatrix(&matrix[i * numOfVariables + i], 1, 1);
		PrintMatrix(eigenVector, 1, numOfEquations);
	}

	file.close();
	return 0;
}




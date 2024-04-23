/**
 * PCC104 - Projeto e Análise de Algoritmos
 * Departamento de Computação - Universidade Federal de Ouro Preto - MG
 * Professor: Pedro Silva
 * Aluno: Fernando dos Santos Alves Fernandes
 * Trabalho Prático I (TP I) - Tarefa 1: Cálculo da multiplicação de duas matrizes
 **/
 
#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>

using namespace std;

class Matrix {
	public:
		Matrix(int size);
		Matrix(vector<vector<int>> matrixElements);
		//~Matrix();
		void setElement(int i, int j, int value);
		void printMatrix();
		void printMatrix(int n);
		void printMatrix(const vector<vector<int>>& matrix);
		vector<vector<int>> squareMatrixMultiply(Matrix B);
		vector<vector<int>> squareMatrixMultiply(const vector<vector<int>>& A, const vector<vector<int>>& B);
		vector<vector<int>> strassen(Matrix B);
		vector<vector<int>> strassen(const vector<vector<int>>& A, const vector<vector<int>>& B);
	private:
		vector<vector<int>> elements;
		vector<vector<int>> *ptE;
		vector<vector<int>> addMatrices(Matrix B);
		vector<vector<int>> addMatrices(const vector<vector<int>>& A, const vector<vector<int>>& B);
		vector<vector<int>> subtractMatrices(Matrix B);
		vector<vector<int>> subtractMatrices(const vector<vector<int>>& A, const vector<vector<int>>& B);	
		
};

// Método construtor.
Matrix::Matrix(int size){
	// Se a ordem das matrizes não for uma potência de 2, preencha com zeros
	int nextPowerOfTwo = 1;
	while (nextPowerOfTwo < size) {
		nextPowerOfTwo *= 2;
	}

	ptE = new vector<vector<int>>(nextPowerOfTwo, vector<int>(nextPowerOfTwo));
	elements = *ptE;
}

// Método construtor.
Matrix::Matrix(vector<vector<int>> matrixElements){
	elements = matrixElements;
}

void Matrix::setElement(int i, int j, int value){
	elements[i][j] = value;
}

// Método para imprimir uma matriz.
void Matrix::printMatrix() {
    for (const auto& row : elements) {
        for (int val : row) {
            cout << val << " ";
        }
        cout << endl;
    }
}

// Método para imprimir uma matriz.
void Matrix::printMatrix(int n) {
    for (int i = 0; i < n; ++i){
        for (int j = 0; j < n; ++j)
			cout << elements[i][j] << " ";
        cout << endl;
    }
}

// Método para imprimir uma matriz.
void Matrix::printMatrix(const vector<vector<int>>& matrix) {
    for (const auto& row : matrix) {
        for (int val : row) {
            cout << val << " ";
        }
        cout << endl;
    }
}

// Método para adicionar duas matrizes.
vector<vector<int>> Matrix::addMatrices(Matrix B) {
    int n = elements.size();
    vector<vector<int>> C(n, vector<int>(n, 0));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            C[i][j] = elements[i][j] + B.elements[i][j];
    return C;
}

// Método para adicionar duas matrizes.
vector<vector<int>> Matrix::addMatrices(const vector<vector<int>>& A, const vector<vector<int>>& B) {
    int n = A.size();
    vector<vector<int>> C(n, vector<int>(n, 0));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            C[i][j] = A[i][j] + B[i][j];
    return C;
}

// Método para subtrair duas matrizes.
vector<vector<int>> Matrix::subtractMatrices(Matrix B) {
    int n = elements.size();
    vector<vector<int>> C(n, vector<int>(n, 0));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            C[i][j] = elements[i][j] - B.elements[i][j];
    return C;
}

// Método para subtrair duas matrizes.
vector<vector<int>> Matrix::subtractMatrices(const vector<vector<int>>& A, const vector<vector<int>>& B) {
    int n = A.size();
    vector<vector<int>> C(n, vector<int>(n, 0));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            C[i][j] = A[i][j] - B[i][j];
    return C;
}

// Método para multiplicação de matrizes usando Dividir para Conquistar.
vector<vector<int>> Matrix::squareMatrixMultiply(Matrix B) {
    int n = elements.size();
    vector<vector<int>> C(n, vector<int>(n, 0));
    if (n == 1) {
        C[0][0] = elements[0][0] * B.elements[0][0];
    } else {
        // Dividindo as matrizes em submatrizes
        int newSize = n / 2;
				
		Matrix mA11(newSize);
		Matrix mA12(newSize);
		Matrix mA21(newSize);
		Matrix mA22(newSize);
		Matrix mB11(newSize);
		Matrix mB12(newSize);
		Matrix mB21(newSize);
		Matrix mB22(newSize);
		
        for (int i = 0; i < newSize; ++i) {
            for (int j = 0; j < newSize; ++j) {
                mA11.setElement(i, j, elements[i][j]);
                mA12.setElement(i, j, elements[i][j + newSize]);
                mA21.setElement(i, j, elements[i + newSize][j]);
                mA22.setElement(i, j, elements[i + newSize][j + newSize]);

                mB11.setElement(i, j, B.elements[i][j]);
                mB12.setElement(i, j, B.elements[i][j + newSize]);
                mB21.setElement(i, j, B.elements[i + newSize][j]);
                mB22.setElement(i, j, B.elements[i + newSize][j + newSize]);
            }
        }
        // Calculando os elementos da matriz C
        vector<vector<int>> C11 = Matrix(mA11.squareMatrixMultiply(mB11)).addMatrices(Matrix(mA12.squareMatrixMultiply(mB21)));
        vector<vector<int>> C12 = Matrix(mA11.squareMatrixMultiply(mB12)).addMatrices(Matrix(mA12.squareMatrixMultiply(mB22)));
        vector<vector<int>> C21 = Matrix(mA21.squareMatrixMultiply(mB11)).addMatrices(Matrix(mA22.squareMatrixMultiply(mB21)));
        vector<vector<int>> C22 = Matrix(mA21.squareMatrixMultiply(mB12)).addMatrices(Matrix(mA22.squareMatrixMultiply(mB22)));
		
        // Combinando as submatrizes para obter a matriz C
        for (int i = 0; i < newSize; ++i) {
            for (int j = 0; j < newSize; ++j) {
                C[i][j] = C11[i][j];
                C[i][j + newSize] = C12[i][j];
                C[i + newSize][j] = C21[i][j];
                C[i + newSize][j + newSize] = C22[i][j];
            }
        }
    }
    return C;
}

// Método para multiplicação de matrizes usando Dividir para Conquistar.
vector<vector<int>> Matrix::squareMatrixMultiply(const vector<vector<int>>& A, const vector<vector<int>>& B) {
    int n = A.size();
    vector<vector<int>> C(n, vector<int>(n, 0));
    if (n == 1) {
        C[0][0] = A[0][0] * B[0][0];
    } else {
        // Dividindo as matrizes em submatrizes
        int newSize = n / 2;
        vector<vector<int>> A11(newSize, vector<int>(newSize));
        vector<vector<int>> A12(newSize, vector<int>(newSize));
        vector<vector<int>> A21(newSize, vector<int>(newSize));
        vector<vector<int>> A22(newSize, vector<int>(newSize));
        vector<vector<int>> B11(newSize, vector<int>(newSize));
        vector<vector<int>> B12(newSize, vector<int>(newSize));
        vector<vector<int>> B21(newSize, vector<int>(newSize));
        vector<vector<int>> B22(newSize, vector<int>(newSize));
        for (int i = 0; i < newSize; ++i) {
            for (int j = 0; j < newSize; ++j) {
                A11[i][j] = A[i][j];
                A12[i][j] = A[i][j + newSize];
                A21[i][j] = A[i + newSize][j];
                A22[i][j] = A[i + newSize][j + newSize];

                B11[i][j] = B[i][j];
                B12[i][j] = B[i][j + newSize];
                B21[i][j] = B[i + newSize][j];
                B22[i][j] = B[i + newSize][j + newSize];
            }
        }
		
        // Calculando os elementos da matriz C
        vector<vector<int>> C11 = addMatrices(squareMatrixMultiply(A11, B11), squareMatrixMultiply(A12, B21));
        vector<vector<int>> C12 = addMatrices(squareMatrixMultiply(A11, B12), squareMatrixMultiply(A12, B22));
        vector<vector<int>> C21 = addMatrices(squareMatrixMultiply(A21, B11), squareMatrixMultiply(A22, B21));
        vector<vector<int>> C22 = addMatrices(squareMatrixMultiply(A21, B12), squareMatrixMultiply(A22, B22));
		
        // Combinando as submatrizes para obter a matriz C
        for (int i = 0; i < newSize; ++i) {
            for (int j = 0; j < newSize; ++j) {
                C[i][j] = C11[i][j];
                C[i][j + newSize] = C12[i][j];
                C[i + newSize][j] = C21[i][j];
                C[i + newSize][j + newSize] = C22[i][j];
            }
        }
    }
    return C;
}

// Método para multiplicar matrizes usando o algoritmo de Strassen.
vector<vector<int>> Matrix::strassen(Matrix B) {
    int n = elements.size();
    vector<vector<int>> C(n, vector<int>(n, 0));
    if (n == 1) {
        C[0][0] = elements[0][0] * B.elements[0][0];
    } else {
        // Dividindo as matrizes em submatrizes
        int newSize = n / 2;
		
		Matrix mA11(newSize);
		Matrix mA12(newSize);
		Matrix mA21(newSize);
		Matrix mA22(newSize);
		Matrix mB11(newSize);
		Matrix mB12(newSize);
		Matrix mB21(newSize);
		Matrix mB22(newSize);
		
        for (int i = 0; i < newSize; ++i) {
            for (int j = 0; j < newSize; ++j) {
                mA11.setElement(i, j, elements[i][j]);
                mA12.setElement(i, j, elements[i][j + newSize]);
                mA21.setElement(i, j, elements[i + newSize][j]);
                mA22.setElement(i, j, elements[i + newSize][j + newSize]);

                mB11.setElement(i, j, B.elements[i][j]);
                mB12.setElement(i, j, B.elements[i][j + newSize]);
                mB21.setElement(i, j, B.elements[i + newSize][j]);
                mB22.setElement(i, j, B.elements[i + newSize][j + newSize]);
            }
        }

        // Calculando as submatrizes
		vector<vector<int>> P1 = mA11.strassen(Matrix(mB12.subtractMatrices(mB22)));
		vector<vector<int>> P2 = Matrix(mA11.addMatrices(mA12)).strassen(mB22);
		vector<vector<int>> P3 = Matrix(mA21.addMatrices(mA22)).strassen(mB11);
		vector<vector<int>> P4 = mA22.strassen(Matrix(mB21.subtractMatrices(mB11)));      
		vector<vector<int>> P5 = Matrix(mA11.addMatrices(mA22)).strassen(Matrix(mB11.addMatrices(mB22)));
		vector<vector<int>> P6 = Matrix(mA12.subtractMatrices(mA22)).strassen(Matrix(mB21.addMatrices(mB22)));
		vector<vector<int>> P7 = Matrix(mA11.subtractMatrices(mA21)).strassen(Matrix(mB11.addMatrices(mB12)));
		
        // Calculando os elementos da matriz C
        vector<vector<int>> C11 = addMatrices(subtractMatrices(addMatrices(P5, P4), P2), P6);
        vector<vector<int>> C12 = addMatrices(P1, P2);
        vector<vector<int>> C21 = addMatrices(P3, P4);
        vector<vector<int>> C22 = subtractMatrices(subtractMatrices(addMatrices(P5, P1), P3), P7);
		
        // Combinando as submatrizes para obter a matriz C
        for (int i = 0; i < newSize; ++i) {
            for (int j = 0; j < newSize; ++j) {
                C[i][j] = C11[i][j];
                C[i][j + newSize] = C12[i][j];
                C[i + newSize][j] = C21[i][j];
                C[i + newSize][j + newSize] = C22[i][j];
            }
        }
    }
    return C;
}

// Método para multiplicar matrizes usando o algoritmo de Strassen.
vector<vector<int>> Matrix::strassen(const vector<vector<int>>& A, const vector<vector<int>>& B) {
    int n = A.size();
    vector<vector<int>> C(n, vector<int>(n, 0));
    if (n == 1) {
        C[0][0] = A[0][0] * B[0][0];
    } else {
        // Dividindo as matrizes em submatrizes
        int newSize = n / 2;
        vector<vector<int>> A11(newSize, vector<int>(newSize));
        vector<vector<int>> A12(newSize, vector<int>(newSize));
        vector<vector<int>> A21(newSize, vector<int>(newSize));
        vector<vector<int>> A22(newSize, vector<int>(newSize));
        vector<vector<int>> B11(newSize, vector<int>(newSize));
        vector<vector<int>> B12(newSize, vector<int>(newSize));
        vector<vector<int>> B21(newSize, vector<int>(newSize));
        vector<vector<int>> B22(newSize, vector<int>(newSize));
        for (int i = 0; i < newSize; ++i) {
            for (int j = 0; j < newSize; ++j) {
                A11[i][j] = A[i][j];
                A12[i][j] = A[i][j + newSize];
                A21[i][j] = A[i + newSize][j];
                A22[i][j] = A[i + newSize][j + newSize];

                B11[i][j] = B[i][j];
                B12[i][j] = B[i][j + newSize];
                B21[i][j] = B[i + newSize][j];
                B22[i][j] = B[i + newSize][j + newSize];
            }
        }
		
        // Calculando as submatrizes
        vector<vector<int>> P1 = strassen(A11, subtractMatrices(B12, B22));
        vector<vector<int>> P2 = strassen(addMatrices(A11, A12), B22);
        vector<vector<int>> P3 = strassen(addMatrices(A21, A22), B11);
        vector<vector<int>> P4 = strassen(A22, subtractMatrices(B21, B11));
        vector<vector<int>> P5 = strassen(addMatrices(A11, A22), addMatrices(B11, B22));
        vector<vector<int>> P6 = strassen(subtractMatrices(A12, A22), addMatrices(B21, B22));
        vector<vector<int>> P7 = strassen(subtractMatrices(A11, A21), addMatrices(B11, B12));
		
        // Calculando os elementos da matriz C
        vector<vector<int>> C11 = addMatrices(subtractMatrices(addMatrices(P5, P4), P2), P6);
        vector<vector<int>> C12 = addMatrices(P1, P2);
        vector<vector<int>> C21 = addMatrices(P3, P4);
        vector<vector<int>> C22 = subtractMatrices(subtractMatrices(addMatrices(P5, P1), P3), P7);
		
        // Combinando as submatrizes para obter a matriz C
        for (int i = 0; i < newSize; ++i) {
            for (int j = 0; j < newSize; ++j) {
                C[i][j] = C11[i][j];
                C[i][j + newSize] = C12[i][j];
                C[i + newSize][j] = C21[i][j];
                C[i + newSize][j + newSize] = C22[i][j];
            }
        }
    }
    return C;
}

// Função principal.
int main(int argc, char** argv){
    if (argc < 3){
        cout << "Execute: file_to_matrix.exe <arquivo.dat> <algoritmo [m: square-matrix-multiply | s: strassen algorithm]>" << endl;
        return -1;
    }
    //else
        //cout << "Arquivo de entrada: " << argv[1] << endl;

    ifstream inputFile(argv[1]);
    if (!inputFile.is_open()) {
        cerr << "Erro ao abrir o arquivo." << endl;
        return 1;
    }
	
	char algorithm = *argv[2];
	//cout << algorithm << endl;

    int n, value;
    inputFile >> n;

	int isPowerOf2 = n && !(n & (n - 1)); // Verifica se n é potência de 2.
	//cout << isPowerOf2 << endl;

	Matrix A(n);
	Matrix B(n);
	
    // Leitura das matrizes A e B do arquivo
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j){
			inputFile >> value;
			A.setElement(i, j, value);
		}

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j){
			inputFile >> value;
			B.setElement(i, j, value);
		}
	
    inputFile.close();

	if (algorithm == 'm') {
		// Multiplicação das matrizes usando o algoritmo Square Matrix Multiply.
		// Início medição de tempo
		auto begin = chrono::high_resolution_clock::now();
		Matrix C(A.squareMatrixMultiply(B));
		// Fim medição de tempo
		auto end = chrono::high_resolution_clock::now();
		auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin);
		//cout << "\nTempo de execucao [Square Matrix Multiply]: " << elapsed.count() * 1e-9 << endl;
		// Impressão da matriz resultante
		//cout << "\nMatrix C:" << endl;
		if (isPowerOf2)
			C.printMatrix();
		else
			C.printMatrix(n);
	} else {
		// Multiplicação das matrizes usando o algoritmo de Strassen.
		// Início medição de tempo
		auto begin = chrono::high_resolution_clock::now();
		Matrix C(A.strassen(B));
		// Fim medição de tempo
		auto end = chrono::high_resolution_clock::now();
		auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin);
		//cout << "\nTempo de execucao [Algoritmo de Strassen]: " << elapsed.count() * 1e-9 << endl;
		// Impressão da matriz resultante
		//cout << "\nMatrix C:" << endl;
		if (isPowerOf2)
			C.printMatrix();
		else
			C.printMatrix(n);
	}
    return 0;
}

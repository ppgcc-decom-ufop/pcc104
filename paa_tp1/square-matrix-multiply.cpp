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

// Função para adicionar duas matrizes.
vector<vector<int>> addMatrices(const vector<vector<int>>& A, const vector<vector<int>>& B) {
    int n = A.size();
    vector<vector<int>> C(n, vector<int>(n, 0));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            C[i][j] = A[i][j] + B[i][j];
    return C;
}

// Função para subtrair duas matrizes.
vector<vector<int>> subtractMatrices(const vector<vector<int>>& A, const vector<vector<int>>& B) {
    int n = A.size();
    vector<vector<int>> C(n, vector<int>(n, 0));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            C[i][j] = A[i][j] - B[i][j];
    return C;
}

// Função para imprimir uma matriz.
void printMatrix(const vector<vector<int>>& matrix) {
    for (const auto& row : matrix) {
        for (int val : row) {
            cout << val << " ";
        }
        cout << endl;
    }
}

// Função para multiplicar matrizes usando o algoritmo de Strassen.
vector<vector<int>> strassen(const vector<vector<int>>& A, const vector<vector<int>>& B) {
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

// Função para multiplicação de matrizes usando Dividir para Conquistar.
vector<vector<int>> squareMatrixMultiply(const vector<vector<int>>& A, const vector<vector<int>>& B) {
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

// Função principal.
int main(int argc, char** argv){
    if (argc < 2){
		cout << "Execute: file_to_matrix.exe <arquivo.dat>" << endl;
		return -1;
	}
	else
		cout << "Arquivo de entrada: " << argv[1] << endl;

	ifstream inputFile(argv[1]);
    if (!inputFile.is_open()) {
        cerr << "Erro ao abrir o arquivo." << endl;
        return 1;
    }

    int n;
    inputFile >> n;
    vector<vector<int>> A(n, vector<int>(n));
    vector<vector<int>> B(n, vector<int>(n));

    // Leitura das matrizes A e B do arquivo
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            inputFile >> A[i][j];

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            inputFile >> B[i][j];

    inputFile.close();
	
	cout << "\nMatriz A:" << endl;
	printMatrix(A);
	
	cout << "\nMatriz B:" << endl;
	printMatrix(B);

    // Multiplicação das matrizes usando o algoritmo Square Matrix Multiply.
	// Início medição de tempo
	auto begin = chrono::high_resolution_clock::now();
	vector<vector<int>> C1 = squareMatrixMultiply(A, B);
	// Fim medição de tempo
	auto end = chrono::high_resolution_clock::now();
	auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin);
	cout << "\nTempo de execucao [Square Matrix Multiply]: " << elapsed.count() * 1e-9 << endl;
	// Impressão da matriz resultante
	cout << "Matrix C:" << endl;
    printMatrix(C1);

    // Multiplicação das matrizes usando o algoritmo de Strassen
	// Início medição de tempo
	begin = chrono::high_resolution_clock::now();
	vector<vector<int>> C2 = strassen(A, B);
	// Fim medição de tempo
	end = chrono::high_resolution_clock::now();
	elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin);
	cout << "\nTempo de execucao [Algoritmo de Strassen]: " << elapsed.count() * 1e-9 << endl;
    // Impressão da matriz resultante
    cout << "Matrix C:" << endl;
    printMatrix(C2);

    return 0;
}

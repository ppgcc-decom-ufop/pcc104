/**
 * PCC104 - Projeto e Analise de Algoritmos
 * Departamento de Computacao - Universidade Federal de Ouro Preto - MG
 * Professor: Pedro Silva
 * Aluno: Fernando dos Santos Alves Fernandes
 * Trabalho Pratico II (TP II) - Tarefa 1: Problema do Clique Maximo
 **/

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>
#include <bits/stdc++.h> //setprecision

using namespace std;

class Graph {
private:
    int V;
    vector<vector<int> > adj;
    vector<int> maxClique;
    int maxCliqueSize;
	double maxTime;
	long numOfCalls = 0;
	std::chrono::time_point<std::chrono::_V2::system_clock, std::chrono::duration<long long int, std::ratio<1ll, 1000000000ll> > > initTime;
	std::chrono::time_point<std::chrono::_V2::system_clock, std::chrono::duration<long long int, std::ratio<1ll, 1000000000ll> > > endTime;
	std::chrono::duration<long long int, std::ratio<1ll, 1000000000ll> > elapsedTime;

public:	
    Graph(int V) : V(V) {
        adj.resize(V);
    }

    void addEdge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }	
	
    // Bron-Kerbosch algorithm (Backtracking version)
	void bronKerboschBT(vector<int>& R, vector<int>& P, vector<int>& X, Graph& g, vector<int>& maxClique) {
		numOfCalls++;
		endTime = chrono::high_resolution_clock::now();
		auto elapsed = chrono::duration_cast<chrono::nanoseconds>(endTime - initTime);
		
		if ((P.empty() && X.empty()) || (elapsed >= chrono::duration_cast<chrono::nanoseconds>(std::chrono::duration<double>(maxTime)))) {
			elapsedTime = elapsed;
			if (R.size() > maxClique.size()) {
				maxClique = R;
			}
			return;
		}

		vector<int> P_copy = P;

		for (int v : P_copy) {
			vector<int> R_new = R; // R'
			R_new.push_back(v); // R' = R U {v}

			vector<int> P_new; // P'
			vector<int> X_new; // X'

			for (int w : P) { // P' = P intersection N(v)
				if (find(g.adj[v].begin(), g.adj[v].end(), w) != g.adj[v].end()) {
					P_new.push_back(w);
				}
			}

			for (int w : X) { // X' = X intersection N(v)
				if (find(g.adj[v].begin(), g.adj[v].end(), w) != g.adj[v].end()) {
					X_new.push_back(w);
				}
			}

			bronKerboschBT(R_new, P_new, X_new, g, maxClique); // bronKerboschBT(R U {v}, P intersection N(v), X intersection N(v)

			P.erase(find(P.begin(), P.end(), v)); // P = P \ {v}
			X.push_back(v); // X = X U {v}
		}
	}
	
	// Bron-Kerbosch algorithm with pivot and bounding (Branch&Bound version)
	void bronKerboschBB(vector<int>& R, vector<int>& P, vector<int>& X, Graph& g, vector<int>& maxClique, int& maxCliqueSize) {
		numOfCalls++;
		endTime = chrono::high_resolution_clock::now();
		auto elapsed = chrono::duration_cast<chrono::nanoseconds>(endTime - initTime);
		//cout << elapsed.count() * 1e-9 << ", " << maxTime << endl;
		
		if ((P.empty() && X.empty()) || (elapsed >= chrono::duration_cast<chrono::nanoseconds>(std::chrono::duration<double>(maxTime)))) {
			if (R.size() > maxCliqueSize) {
				elapsedTime = elapsed;
				maxCliqueSize = R.size();
				maxClique = R;
			}
			return;
		}

		int pivot = P.empty() ? -1 : P[0]; // Se P for vazio, nao ha como escolher um pivo.
		vector<int> P_aux; // Candidatos que nao sao vizinhos de v.
		for (int v : P) {
			if (pivot == -1 || find(g.adj[pivot].begin(), g.adj[pivot].end(), v) == g.adj[pivot].end()) {
				P_aux.push_back(v); // Vertices que nao sao adjacentes a pivot.
			}
		}

		for (int v : P_aux) {
			vector<int> R_new = R;
			R_new.push_back(v);

			vector<int> P_new;
			vector<int> X_new;

			for (int w : P) {
				if (find(g.adj[v].begin(), g.adj[v].end(), w) != g.adj[v].end()) {
					P_new.push_back(w);
				}
			}

			for (int w : X) {
				if (find(g.adj[v].begin(), g.adj[v].end(), w) != g.adj[v].end()) {
					X_new.push_back(w);
				}
			}

			// So prossegue se o tamanho potencial do clique e maior do que o clique maximo atualmente encontrado.
			if (R_new.size() + P_new.size() > maxCliqueSize) {
				bronKerboschBB(R_new, P_new, X_new, g, maxClique, maxCliqueSize);
			}

			P.erase(find(P.begin(), P.end(), v));
			X.push_back(v);
		}
	}	
	
	void setInitialAndMaxTime(std::chrono::time_point<std::chrono::_V2::system_clock, std::chrono::duration<long long int, std::ratio<1ll, 1000000000ll> > > initial, double time) {
		initTime = initial;
		maxTime = time;
	}
	
	auto getElasedTime(){
		return elapsedTime.count();
	}
	
	void printElapsedTime(){
		cout << "\nTempo de execucao: " << std::fixed << std::setprecision(4) << double (elapsedTime.count() * 1e-9) << " s."<< endl;
	}
	
    void printMaxClique() {
		cout << maxCliqueSize << endl;

        for ( int v : maxClique) {
            cout << v + 1 << " ";
        }
        cout << endl;
    }
	
	long getNumOfCalls(){
		return numOfCalls;
	}
};

int main(int argc, char* argv[]) {
    if (argc != 4) {
        cerr << "Uso: " << argv[0] << " <arquivoProblema> <tempoMaximo> <estrategia>" << endl;
        return 1;
    }

    string filename = argv[1];
    int timetorun = stoi(argv[2]);
	string strategy = argv[3];

    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Erro ao abrir o arquivo " << filename << endl;
        return 1;
    }

    string line;
    unsigned int numVertices = 0;
    unsigned int numEdges = 0;

	string extension = ".grafo";
	string::size_type i = filename.find(extension);
	if (i != string::npos) { // extensao '.grafo'
		getline(file, line);
		stringstream ss(line);
		ss >> numVertices;
		ss >> numEdges;
	}else {	// extensao '.clq'
		// Ler o cabecalho do arquivo e obter o numero de vertices e arestas
		while (getline(file, line)) {
			if (line.empty() || line[0] == 'p') {
				stringstream ss(line);
				string token;
				ss >> token; // descartar 'p'
				ss >> token; // esperado 'col'/'edge'
				ss >> numVertices;
				ss >> numEdges;
				break;
			}
		}
	}
	
    Graph graph(numVertices);

    // Ler as arestas e adicionar ao grafo
	if (i != string::npos) { // extensao '.grafo'
		while (getline(file, line)) {
			if (line.empty()) continue;
			unsigned int from, to;
			stringstream ss(line);
			ss >> from >> to;
			graph.addEdge(from, to); // Os vertices no arquivo comecam de 0.
		}
	} else{
		while (getline(file, line)) {
			if (line.empty() || line[0] != 'e') continue;
			stringstream ss(line);
			string token;
			ss >> token; // descartar 'e'
			unsigned int from, to;
			ss >> from >> to;
			graph.addEdge(from - 1, to - 1); // Os vertices no arquivo comecam de 1.
		}
	}
    file.close();

    vector<int> R, P, X, maxClique;
	int maxCliqueSize = 0;
    for (int i = 0; i < numVertices; ++i) {
        P.push_back(i);
    }
	
	auto initTime = chrono::high_resolution_clock::now();
	graph.setInitialAndMaxTime(initTime, timetorun);
	if (strategy == "bp") // Backtracking ou Backpropagation
		graph.bronKerboschBT(R, P, X, graph, maxClique);
	else if (strategy == "bb") // Branch and Bound
		graph.bronKerboschBB(R, P, X, graph, maxClique, maxCliqueSize);
    //graph.printMaxClique();
	
	cout << "Tempo de execucao: " << graph.getElasedTime() * 1e-9 << endl;
	cout << "Numero de chamadas recursivas: " << graph.getNumOfCalls() << endl;
	
	cout << maxClique.size() << endl;
    for (int v : maxClique) {
        cout << v + 1 << " "; // Convert back to 1-based index
    }
    cout << endl;
	//graph.printElapsedTime();

    return 0;
}
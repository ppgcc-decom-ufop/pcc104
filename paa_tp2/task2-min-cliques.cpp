/**
 * PCC104 - Projeto e Analise de Algoritmos
 * Departamento de Computacao - Universidade Federal de Ouro Preto - MG
 * Professor: Pedro Silva
 * Aluno: Fernando dos Santos Alves Fernandes
 * Trabalho Pratico II (TP II) - Tarefa 2: Problema do Minimo Numero de Cliques
 **/

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <limits>
#include <chrono>
#include <algorithm>
#include <bits/stdc++.h> //setprecision
#include <unordered_set>

using namespace std;

class Graph {
private:
    int V;
    vector<vector<int> > adj;
    vector<vector<int> > bestCover; // Solucao
    int bestCoverSize;
    double maxTime;
	long numOfCalls = 0;
    chrono::time_point<chrono::system_clock> initTime;

public:
    Graph(int V) : V(V), bestCoverSize(numeric_limits<int>::max()) {
        adj.resize(V);
    }

    void addEdge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    void setInitialAndMaxTime(chrono::time_point<chrono::system_clock> initial, double time) {
        initTime = initial;
        maxTime = time;
    }

	// Algoritmo de Bron-Kerbosch para encontrar todos os cliques maximais (versao backtracking).
    void findSetOfCliques(vector<int>& R, vector<int>& P, vector<int>& X, vector<vector<int> >& allCliques) {
		numOfCalls++;
        if (P.empty() && X.empty()) {
            allCliques.push_back(R);
            return;
        }

        while (!P.empty()) {
            int v = P.back();
            P.pop_back();
            vector<int> R_new = R;
            R_new.push_back(v);

            vector<int> P_new;
            vector<int> X_new;
            for (int w : P) {
                if (find(adj[v].begin(), adj[v].end(), w) != adj[v].end()) {
                    P_new.push_back(w);
                }
            }
            for (int w : X) {
                if (find(adj[v].begin(), adj[v].end(), w) != adj[v].end()) {
                    X_new.push_back(w);
                }
            }

            findSetOfCliques(R_new, P_new, X_new, allCliques);
            X.push_back(v);
			
			// Limite de tempo de execucao
            if (chrono::duration_cast<chrono::seconds>(chrono::system_clock::now() - initTime).count() > maxTime) {
                return;
            }
        }
    }

    vector<vector<int> > cliqueCover() {
        vector<vector<int> > allCliques;
        vector<int> R, P(V), X;
        iota(P.begin(), P.end(), 0); // Inicializa P (conjunto de vertices candidados) com todos os vertices.

        findSetOfCliques(R, P, X, allCliques);
        return allCliques;
    }

	bool compareSize(vector<int> clique1, vector<int> clique2) {
		return (clique1.size() > clique2.size());
	}

	void printSolution() {
        // Ordena os cliques pelo tamanho, do maior para o menor.
        sort(bestCover.begin(), bestCover.end(), [](const vector<int>& a, const vector<int>& b) {
            return a.size() > b.size();
        });

        cout << bestCover.size() << endl;
        for (const auto& clique : bestCover) {
            for (int v : clique) {
                cout << v + 1 << " ";
            }
            cout << "\n";
        }
    }

    void minimumCliqueCover() {
        auto allCliques = cliqueCover();
        vector<int> covered(V, 0);
        vector<vector<int> > min_cover;

        for (const auto& clique : allCliques) { // Para cada clique, checa se cobre algum vertice ainda nao coberto.
            bool is_covered = false;
            for (int v : clique) {  
                if (!covered[v]) {
                    is_covered = true;
                    break;
                }
            }
            if (is_covered) { // Se pelo menos um vertice ainda nao foi coberto, todos os vertices do clique passam a ser cobertos.
                min_cover.push_back(clique); // O clique e incluido na cobertura.
                for (int v : clique) {
                    covered[v] = 1;
                }
            }
        }

        bestCover = min_cover;
    }
	
	long getNumOfCalls(){
		return numOfCalls;
	}
};

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Uso: " << argv[0] << " <arquivoProblema> <tempoMaximo>" << endl;
        return 1;
    }

    string filename = argv[1];
    double timetorun = stod(argv[2]);

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
    if (i != string::npos) {
        getline(file, line);
        stringstream ss(line);
        ss >> numVertices;
        ss >> numEdges;
    } else {
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

    if (i != string::npos) {
        while (getline(file, line)) {
            if (line.empty()) continue;
            unsigned int from, to;
            stringstream ss(line);
            ss >> from >> to;
            graph.addEdge(from, to);
        }
    } else {
        while (getline(file, line)) {
            if (line.empty() || line[0] != 'e') continue;
            stringstream ss(line);
            string token;
            ss >> token; // descartar 'e'
            unsigned int from, to;
            ss >> from >> to;
            graph.addEdge(from - 1, to - 1);
        }
    }
    file.close();

    auto startTime = chrono::system_clock::now();
    graph.setInitialAndMaxTime(startTime, timetorun);

    graph.minimumCliqueCover();
	//cout << "Numero de chamadas recursivas: " << graph.getNumOfCalls() << endl;
    graph.printSolution();

    return 0;
}

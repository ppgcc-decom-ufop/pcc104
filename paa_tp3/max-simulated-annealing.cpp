/**
 * PCC104 - Projeto e Analise de Algoritmos
 * Departamento de Computacao - Universidade Federal de Ouro Preto - MG
 * Professor: Pedro Silva
 * Aluno: Fernando dos Santos Alves Fernandes
 * Trabalho Pratico III (TP III) - Tarefa 1: Problema do Clique Maximo
 * Algoritmo: Simulated Annealing
 **/

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>
#include <bits/stdc++.h> //setprecision
#include <ctime>
#include <cstdlib>
#include <ctime>
#include <cstdlib>
#include <cmath>

using namespace std;

class Graph {
public:
    int V;
    std::vector<std::vector<int>> adj;

    Graph(int V) : V(V), adj(V) {}

    void addEdge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    bool isClique(const std::vector<int>& vertices) {
        for (size_t i = 0; i < vertices.size(); ++i) {
            for (size_t j = i + 1; j < vertices.size(); ++j) {
                if (std::find(adj[vertices[i]].begin(), adj[vertices[i]].end(), vertices[j]) == adj[vertices[i]].end()) {
                    return false;
                }
            }
        }
        return true;
    }

    std::vector<int> getRandomClique() {
        std::vector<int> clique;
        int start = rand() % V;
        clique.push_back(start);
        return clique;
    }

    std::vector<int> getNeighbors(const std::vector<int>& clique) {
        std::vector<int> neighbors;
        for (int v : clique) {
            for (int u : adj[v]) {
                if (std::find(clique.begin(), clique.end(), u) == clique.end()) {
                    neighbors.push_back(u);
                }
            }
        }
        return neighbors;
    }

    double probability(double delta, double temperature) {
        return exp(delta / temperature);
    }

    std::vector<int> simulatedAnnealing(double initialTemperature, double coolingRate) {
        std::vector<int> current = getRandomClique();
        std::vector<int> best = current;
        double temperature = initialTemperature;

        while (temperature > 1) {
            std::vector<int> neighbors = getNeighbors(current);
            if (neighbors.empty()) {
                break;
            }
            int idx = rand() % neighbors.size();
            std::vector<int> newClique = current;
            newClique.push_back(neighbors[idx]);

            if (isClique(newClique)) {
                if (newClique.size() > current.size() || probability(newClique.size() - current.size(), temperature) > (double)rand() / RAND_MAX) {
                    current = newClique;
                }
                if (current.size() > best.size()) {
                    best = current;
                }
            }
            temperature *= coolingRate;
        }
        return best;
    }
};

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Uso: " << argv[0] << " <arquivoProblema> <tempoMaximo>" << endl;
        return 1;
    }

    string filename = argv[1];
    int timetorun = stoi(argv[2]);
    //string strategy = argv[3];

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
    }else {    // extensao '.clq'
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
    srand(time(0));
	
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

    std::vector<int> clique = graph.simulatedAnnealing(1000, 0.95);
    //std::cout << "Clique máximo encontrado: ";
	cout << clique.size() << endl;
    for (int v : clique) {
        std::cout << v + 1 << " ";
    }
    std::cout << std::endl;
    return 0;
}

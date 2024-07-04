/**
 * PCC104 - Projeto e Analise de Algoritmos
 * Departamento de Computacao - Universidade Federal de Ouro Preto - MG
 * Professor: Pedro Silva
 * Aluno: Fernando dos Santos Alves Fernandes
 * Trabalho Pratico III (TP III) - Tarefa 1: Problema do Clique Maximo
 * Algoritmo: Hill Climbing
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

    std::vector<int> hillClimbing() {
        std::vector<int> current = getRandomClique();
        while (true) {
            std::vector<int> neighbors = getNeighbors(current);
            std::vector<int> bestNeighbor = current;
            for (int neighbor : neighbors) {
                std::vector<int> newClique = current;
                newClique.push_back(neighbor);
                if (isClique(newClique) && newClique.size() > bestNeighbor.size()) {
                    bestNeighbor = newClique;
                }
            }
            if (bestNeighbor.size() == current.size()) {
                break;
            }
            current = bestNeighbor;
        }
        return current;
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
	
/*     Graph g(7);
    g.addEdge(0, 1);
    g.addEdge(0, 2);
    g.addEdge(1, 2);
    g.addEdge(1, 3);
    g.addEdge(3, 4);
    g.addEdge(3, 5);
	g.addEdge(3, 6);
	g.addEdge(4, 5);
    g.addEdge(4, 6);
    g.addEdge(5, 6); */

    std::vector<int> clique = graph.hillClimbing();
    //std::cout << "Clique máximo encontrado: ";
	cout << clique.size() << endl;
    for (int v : clique) {
        std::cout << v + 1 << " ";
    }
    std::cout << std::endl;
    return 0;
}

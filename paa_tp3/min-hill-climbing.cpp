/**
 * PCC104 - Projeto e Analise de Algoritmos
 * Departamento de Computacao - Universidade Federal de Ouro Preto - MG
 * Professor: Pedro Silva
 * Aluno: Fernando dos Santos Alves Fernandes
 * Trabalho Pratico III (TP III) - Tarefa 2: Problema do Minimo Numero de Cliques
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
#include <unordered_set>
#include <unordered_map>

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

    std::vector<int> getRandomClique(const std::unordered_set<int>& coveredVertices) {
        std::vector<int> clique;
        int start;
        do {
            start = rand() % V;
        } while (coveredVertices.find(start) != coveredVertices.end());
        clique.push_back(start);
        return clique;
    }

    std::vector<int> getNeighbors(const std::vector<int>& clique, const std::unordered_set<int>& coveredVertices) {
        std::vector<int> neighbors;
        for (int v : clique) {
            for (int u : adj[v]) {
                if (std::find(clique.begin(), clique.end(), u) == clique.end() && coveredVertices.find(u) == coveredVertices.end()) {
                    neighbors.push_back(u);
                }
            }
        }
        return neighbors;
    }

    std::vector<int> hillClimbing(const std::unordered_set<int>& coveredVertices) {
        std::vector<int> current = getRandomClique(coveredVertices);
        while (true) {
            std::vector<int> neighbors = getNeighbors(current, coveredVertices);
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

    bool isSubset(const std::vector<int>& clique, const std::vector<std::vector<int>>& cliques) {
        for (const auto& existingClique : cliques) {
            std::unordered_set<int> existingSet(existingClique.begin(), existingClique.end());
            bool isSubset = true;
            for (int v : clique) {
                if (existingSet.find(v) == existingSet.end()) {
                    isSubset = false;
                    break;
                }
            }
            if (isSubset) {
                return true;
            }
        }
        return false;
    }

    std::vector<std::vector<int>> vertexCoverByCliques(int timeLimit) {
        auto start = chrono::steady_clock::now();
        std::vector<std::vector<int>> cliques;
        std::unordered_set<int> coveredVertices;

        while (coveredVertices.size() < V) {
            std::vector<int> clique = hillClimbing(coveredVertices);
            if (!isSubset(clique, cliques)) {
                cliques.push_back(clique);
                for (int v : clique) {
                    coveredVertices.insert(v);
                }
            }

            auto end = chrono::steady_clock::now();
            auto duration = chrono::duration_cast<chrono::seconds>(end - start).count();
            if (duration >= timeLimit) {
                break;
            }
        }
		
		// Ordenar cliques do maior para o menor
        std::sort(cliques.begin(), cliques.end(), [](const std::vector<int>& a, const std::vector<int>& b) {
            return a.size() > b.size();
        });

        return cliques;
    }
};

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Uso: " << argv[0] << " <arquivoProblema> <tempoMaximo>" << endl;
        return 1;
    }

    string filename = argv[1];
    int timetorun = stoi(argv[2]);

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
    } else { // extensao '.clq'
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
    } else {
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

    std::vector<std::vector<int>> cliques = graph.vertexCoverByCliques(timetorun);
    cout << cliques.size() << endl;
    for (const auto& clique : cliques) {
        for (int v : clique) {
            cout << v + 1 << " ";
        }
        cout << endl;
    }

    return 0;
}

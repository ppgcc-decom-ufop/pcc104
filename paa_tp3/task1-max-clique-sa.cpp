#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>
#include <random>
#include <cmath>
#include <bits/stdc++.h> //setprecision

using namespace std;

class Graph {
private:
    int V;
    vector<vector<int>> adj;
    vector<int> maxClique;
    int maxCliqueSize;
    double maxTime;
    long numOfCalls = 0;
    chrono::time_point<chrono::high_resolution_clock> initTime;
    chrono::time_point<chrono::high_resolution_clock> endTime;
    chrono::duration<double> elapsedTime;

public:
    Graph(int V) : V(V) {
        adj.resize(V);
    }

    void addEdge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    bool isClique(const vector<int>& subset) {
        for (size_t i = 0; i < subset.size(); ++i) {
            for (size_t j = i + 1; j < subset.size(); ++j) {
                int u = subset[i];
                int v = subset[j];
                if (find(adj[u].begin(), adj[u].end(), v) == adj[u].end()) {
                    return false;
                }
            }
        }
        return true;
    }

    vector<int> getRandomClique() {
        vector<int> clique;
        int vertex = rand() % V;
        clique.push_back(vertex);
        for (int i = 0; i < adj[vertex].size(); ++i) {
            clique.push_back(adj[vertex][i]);
        }
        return clique;
    }

    vector<int> perturbSolution(const vector<int>& currentClique) {
        vector<int> newClique = currentClique;
        int action = rand() % 2;
        if (action == 0 && !newClique.empty()) { 
            // Remover um vértice aleatório
            int index = rand() % newClique.size();
            newClique.erase(newClique.begin() + index);
        } else { 
            // Adicionar um vértice aleatório
            int vertex = rand() % V;
            newClique.push_back(vertex);
            // Garantir que o novo conjunto ainda seja um clique
            if (!isClique(newClique)) {
                newClique.pop_back();
            }
        }
        return newClique;
    }

    vector<int> simulatedAnnealing(double initialTemp, double finalTemp, double alpha, int maxIterations) {
        vector<int> currentClique = getRandomClique();
        vector<int> bestClique = currentClique;
        double T = initialTemp;
        int iteration = 0;

        while (T > finalTemp && iteration < maxIterations) {
            vector<int> newClique = perturbSolution(currentClique);

            int currentCost = currentClique.size();
            int newCost = newClique.size();

            if (newCost > currentCost) {
                currentClique = newClique;
                if (newCost > bestClique.size()) {
                    bestClique = newClique;
                }
            } else {
                double acceptanceProbability = exp((newCost - currentCost) / T);
                if (((double) rand() / RAND_MAX) < acceptanceProbability) {
                    currentClique = newClique;
                }
            }

            T *= alpha;
            iteration++;
        }

        return bestClique;
    }

    void setInitialAndMaxTime(chrono::time_point<chrono::high_resolution_clock> initial, double time) {
        initTime = initial;
        maxTime = time;
    }

    auto getElapsedTime() {
        return elapsedTime.count();
    }

    void printElapsedTime() {
        cout << "\nTempo de execucao: " << fixed << setprecision(4) << double(elapsedTime.count()) << " s." << endl;
    }

    long getNumOfCalls() {
        return numOfCalls;
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
                ss >> token;
                ss >> token;
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
            ss >> token;
            unsigned int from, to;
            ss >> from >> to;
            graph.addEdge(from - 1, to - 1);
        }
    }
    file.close();

    vector<int> maxClique;
	double initialTemp = 1000.0;   // Temperatura inicial
	double finalTemp = 0.1;        // Temperatura final
	double alpha = 0.99;           // Taxa de resfriamento
	int maxIterations = 10000;     // Número máximo de iterações

	auto initTime = chrono::high_resolution_clock::now();
	graph.setInitialAndMaxTime(initTime, timetorun);
	maxClique = graph.simulatedAnnealing(initialTemp, finalTemp, alpha, maxIterations);

    //cout << "Tempo de execucao: " << graph.getElapsedTime() * 1e-9 << endl;
    //cout << "Numero de chamadas recursivas: " << graph.getNumOfCalls() << endl;
    //cout << "Tamanho do clique maximo: " << maxClique.size() << endl;
	cout << maxClique.size() << endl;
    for (int v : maxClique) {
        cout << v + 1 << " "; // Convert back to 1-based index
    }
    cout << endl;

    return 0;
}

/**
 * PCC104 - Projeto e Analise de Algoritmos
 * Departamento de Computacao - Universidade Federal de Ouro Preto - MG
 * Professor: Pedro Silva
 * Aluno: Fernando dos Santos Alves Fernandes
 * Trabalho Pratico III (TP III) - Tarefa 1: Problema do Clique Maximo
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

    bool viavel(vector<int>& clique) {
        if (clique.size() < 1)
            return false;
        else {
            if (clique.size() == 1)
                return true;
            for (int v: clique){
                for (int w: clique) {
                    if (find(adj[v].begin(), adj[v].end(), w) == adj[v].end())
                        return false;
                }
            }
        }
    }

    bool isClique(const vector<int>& subset) {
        for (size_t i = 0; i < subset.size(); ++i) {
            for (size_t j = i + 1; j < subset.size(); ++j) {
                int u = subset[i];
                int v = subset[j];
                if (find(adj[u].begin(), adj[u].end(), v) == adj[u].end()) {
                    return false; // Se não houver aresta, não é um clique
                }
            }
        }
        return true; // Se todas as arestas existirem, é um clique
    }

    void printMaxClique() {
        cout << maxCliqueSize << endl;

        for ( int v : maxClique) {
            cout << v + 1 << " ";
        }
        cout << endl;
    }

    void printSet(vector<int>& setToPrint) {
        cout << "S: " << setToPrint.size() << " -> ";

        for ( int v : setToPrint) {
            cout << v + 1 << " ";
        }
        cout << endl;
    }

    vector<vector<int>> getSortedCandidates() {
        // Vetor para armazenar os vértices e seus graus
        vector<std::vector<int>> sortedCandidates;

        // Calcular os graus dos vértices
        for (int i = 0; i < adj.size(); ++i) {
            int grau = adj[i].size();
            sortedCandidates.push_back({i, grau});
        }

        // Ordenar os vértices por grau em ordem crescente
        sort(sortedCandidates.begin(), sortedCandidates.end(),
                  [](const vector<int>& a, const vector<int>& b) {
                      return a[1] < b[1];
                  });

        return sortedCandidates;
    }

    void maiorClique(vector<int>& S, vector<int>& C, vector<int>& maxClique) {
        numOfCalls++;
        endTime = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::nanoseconds>(endTime - initTime);
        
        vector<vector<int> > sortedCandidates = getSortedCandidates();
        
        // Ordena os vertices candidatos pelo grau (numero de vertices adjacentes).
        sort(sortedCandidates.begin(), sortedCandidates.end(), [](const vector<int>& a, const vector<int>& b) {
            return (a[1] < b[1]);
        });        
        
        while ((!sortedCandidates.empty()) && (elapsed < chrono::duration_cast<chrono::nanoseconds>(std::chrono::duration<double>(maxTime)))) {
            elapsedTime = elapsed;

            int x = sortedCandidates.back()[0]; // x = seleciona(C);
            
            sortedCandidates.pop_back(); // C = C \ {x};

            S.push_back(x);
            //bool feasible = isClique(S);
            //if (!feasible) { // se nao eh viavel(S U {x}) entao
            if (!isClique(S)) {
                S.pop_back(); // S = S \ {x};
            }
            if (S.size() > maxClique.size()) {
                maxClique = S;
            }
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

    vector<int> S, C, X, maxClique;
    int maxCliqueSize = 0;
    for (int i = 0; i < numVertices; ++i) {
        C.push_back(i);
    }
    
    auto initTime = chrono::high_resolution_clock::now();
    graph.setInitialAndMaxTime(initTime, timetorun);
    graph.maiorClique(S, C, maxClique);
    //graph.printMaxClique();

    //cout << "Tempo de execucao: " << graph.getElasedTime() * 1e-9 << endl;
    //cout << "Numero de chamadas recursivas: " << graph.getNumOfCalls() << endl;
    //cout << "Tamanho do clique maximo: " << maxClique.size() << endl;
    cout << maxClique.size() << endl;
    for (int v : maxClique) {
        cout << v + 1 << " "; // Convert back to 1-based index
    }
    cout << endl;
    //graph.printElapsedTime();

    return 0;
}
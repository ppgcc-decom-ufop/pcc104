/**
 * PCC104 - Projeto e Análise de Algoritmos
 * Departamento de Computação - Universidade Federal de Ouro Preto - MG
 * Professor: Pedro Silva
 * Aluno: Fernando dos Santos Alves Fernandes
 * Trabalho Prático I (TP I) - Tarefa 2: Algoritmo de Dijkstra
 **/

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <limits>
#include <queue>
#include <algorithm>
//#include <chrono>

using namespace std;

// Estrutura para representar uma aresta
struct Edge {
    unsigned int to;
    unsigned short int weight;
    Edge(unsigned int t, unsigned short int w) : to(t), weight(w) {}
};

// Estrutura para representar um grafo
class Graph {
private:
    unsigned int numVertices;
    vector<vector<Edge> > adjacencyList;

public:
    Graph(unsigned int n) : numVertices(n), adjacencyList(n) {}

    // Adiciona uma aresta ao grafo
    void addEdge(unsigned int from, unsigned int to, unsigned short int weight) {
        adjacencyList[from].emplace_back(to, weight);
    }
    
    void printGraph() {
        for (unsigned int v = 0; v < numVertices; v++) {
            cout << "[" << (v + 1) << "]: {";
            for (const Edge& edge : adjacencyList[v]) {
                cout << "(" << (edge.to + 1) << ", " << edge.weight << ")";
            }
            cout << "}" << endl;
        }
    }

    // Função para encontrar o caminho mínimo usando o algoritmo de Dijkstra
    pair<vector<unsigned int>, vector<unsigned int> > dijkstra(unsigned int start) {
        vector<unsigned int> distance(numVertices, numeric_limits<int>::max());
        vector<unsigned int> predecessor(numVertices, -1);
        
        distance[start] = 0;

        priority_queue<pair<int, int>, vector<pair<int, int> >, greater<pair<int, int> >> pq;
        pq.push({0, start});
				
        while (!pq.empty()) {
            unsigned int u = pq.top().second;
            unsigned int dist = pq.top().first;
            pq.pop();
			
            if (dist > distance[u]) continue;


            //cout << "u = " << u << endl;
            for (const Edge& edge : adjacencyList[u]) {
                unsigned int v = edge.to;
                unsigned short int weight = edge.weight;
                //cout << "v = " << v << endl;
                if (distance[v] > distance[u] + weight) {
                    distance[v] = distance[u] + weight;
                    predecessor[v] = u;
                    pq.push({distance[v], v});
                }
            }
        }

        return {distance, predecessor};
    }
    
    vector<unsigned int> reconstructPath(unsigned int source, unsigned int destination, const vector<unsigned int>& predecessor) {
        vector<unsigned int> path;
		//path.reserve(numVertices);
        unsigned int current = destination;
        
		while (current != source) {
			path.push_back(current + 1);
			current = predecessor[current];
			//cout << "current = " << current << endl;
		}
		path.push_back(source + 1);
		//for (int i = 0; i < predecessor.size(); i++)
		//    cout << " " << predecessor[i] << " ";
		//cout << endl;
		reverse(path.begin(), path.end());
		
        return path;
    }
};

int main(int argc, char* argv[]) {
    if (argc != 4) {
        cerr << "Uso: " << argv[0] << " <arquivoProblema> <nrExecuções> <flagDepuração>" << endl;
        return 1;
    }

    string filename = argv[1];
    unsigned short int numExecutions = stoi(argv[2]);
    unsigned short int debugFlag = stoi(argv[3]);

    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Erro ao abrir o arquivo " << filename << endl;
        return 1;
    }

    string line;
    unsigned int numVertices = 0;
    unsigned int numEdges = 0;

    // Ler o cabeçalho do arquivo e obter o número de vértices e arestas
    while (getline(file, line)) {
        if (line.empty() || line[0] == 'p') {
            stringstream ss(line);
            string token;
            ss >> token; // descartar 'p'
            ss >> token; // esperado 'sp'
            ss >> numVertices;
            ss >> numEdges;
            break;
        }
    }

    Graph graph(numVertices);

    // Ler as arestas e adicionar ao grafo
    while (getline(file, line)) {
        if (line.empty() || line[0] != 'a') continue;
        stringstream ss(line);
        string token;
        ss >> token; // descartar 'a'
        unsigned int from, to, weight;
        ss >> from >> to >> weight;
        graph.addEdge(from - 1, to - 1, weight); // Os vértices no arquivo começam de 1
    }
    
    file.close();
    
    //graph.printGraph();
    
    unsigned int startVertex;
    if (debugFlag == 1){
        startVertex = 0;
        numExecutions = 1;
    } else {
        startVertex = rand() % numVertices;
    }
    // Executar Dijkstra para cada vértice
    for (unsigned short int i = 0; i < numExecutions; ++i) {
        // Imprimir os caminhos mínimos em arquivo, se a flag de depuração estiver ativada.
        if (debugFlag) {
	    //auto begin = chrono::high_resolution_clock::now();
	    pair<vector<unsigned int>, vector<unsigned int> > result = graph.dijkstra(startVertex);
	    //auto end = chrono::high_resolution_clock::now();
	    //auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin);
	    //cout << "\nTempo de execucao [graph.dijkstra]: " << elapsed.count() * 1e-9 << endl;
            // Reconstrói o caminho mínimo
            //int targetVertex;// = stoi(argv[4]);
             string s = ".gr";
            string::size_type i = filename.find(s);
            if (i != string::npos){
                filename.erase(i, s.length());
                filename = "spaths"+filename+".txt";
            }
            
            ofstream spath(filename);
            if (!spath.is_open()) {
                cerr << "Erro ao abrir o arquivo " << filename << endl;
                return 1;
            } else {
                for (unsigned int targetVertex = 2; targetVertex <= numVertices; targetVertex++) {
                    vector<unsigned int> path = graph.reconstructPath(startVertex, targetVertex - 1, result.second);
                    // Exibe o caminho mínimo no formato especificado
                    spath << "[" << startVertex + 1 << "," << targetVertex << "](" << result.first[targetVertex - 1] << ") ";
                    for (unsigned int i = 0; i < path.size(); ++i) {
                        spath << path[i] << " ";
                    }
                    spath << endl;
                }
            }
            spath.close();
        } else {        
            cout << "\nstartVertex = " << (startVertex + 1) << endl;
	    //auto begin = chrono::high_resolution_clock::now();
	    pair<vector<unsigned int>, vector<unsigned int> > result = graph.dijkstra(startVertex);
	    //auto end = chrono::high_resolution_clock::now();
	    //auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin);			
	    //cout << "\nTempo de execucao [graph.dijkstra]: " << elapsed.count() * 1e-9 << endl;
            for (unsigned int targetVertex = 1; targetVertex <= numVertices; targetVertex++) {
                vector<unsigned int> path = graph.reconstructPath(startVertex, targetVertex - 1, result.second);
                // Exibe o caminho mínimo no formato especificado
                cout << "\n[" << startVertex + 1 << "," << targetVertex << "](" << result.first[targetVertex - 1] << ") ";
                for (unsigned int i = 0; i < path.size(); ++i) {
                    cout << path[i] << " ";
                }
            }                        
            cout << endl;
            startVertex = rand() % numVertices;
        }
    }
    return 0;
}

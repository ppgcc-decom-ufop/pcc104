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

using namespace std;

// Estrutura para representar uma aresta
struct Edge {
    int to;
    int weight;
    Edge(int t, int w) : to(t), weight(w) {}
};

// Estrutura para representar um grafo
class Graph {
private:
    int numVertices;
    vector<vector<Edge>> adjacencyList;

public:
    Graph(int n) : numVertices(n), adjacencyList(n) {}

    // Adiciona uma aresta ao grafo
    void addEdge(int from, int to, int weight) {
        adjacencyList[from].emplace_back(to, weight);
    }
	
	void printGraph() {
		for (int v = 0; v < numVertices; v++) {
			cout << "[" << (v + 1) << "]: {";
			for (const Edge& edge : adjacencyList[v]) {
				cout << "(" << (edge.to + 1) << ", " << edge.weight << ")";
			}
			cout << "}" << endl;
		}
	}

    // Função para encontrar o caminho mínimo usando o algoritmo de Dijkstra
    pair<vector<int>, vector<int>> dijkstra(int start) {
        vector<int> distance(numVertices, numeric_limits<int>::max());
		vector<int> predecessor(numVertices, -1);
		
        distance[start] = 0;

        priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
        pq.push({0, start});

        while (!pq.empty()) {
            int u = pq.top().second;
            int dist = pq.top().first;
            pq.pop();
			
            if (dist > distance[u]) continue;

			//cout << "u = " << u << endl;
            for (const Edge& edge : adjacencyList[u]) {
                int v = edge.to;
                int weight = edge.weight;
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
	
	vector<int> reconstructPath(int source, int destination, const vector<int>& predecessor) {
		vector<int> path;
		int current = destination;
		
		while (current != source) {
			path.push_back(current + 1);
			current = predecessor[current];
			//cout << "current = " << current << endl;
		}
		path.push_back(source + 1);
		//for (int i = 0; i < predecessor.size(); i++)
		//	cout << " " << predecessor[i] << " ";
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
    int numExecutions = stoi(argv[2]);
    int debugFlag = stoi(argv[3]);

    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Erro ao abrir o arquivo " << filename << endl;
        return 1;
    }

    string line;
    int numVertices = 0;
    int numEdges = 0;

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
        int from, to, weight;
        ss >> from >> to >> weight;
        graph.addEdge(from - 1, to - 1, weight); // Os vértices no arquivo começam de 1
    }
	
	file.close();
	
	//graph.printGraph();
	
 
	int startVertex;
	if (debugFlag == 1){
		startVertex = 0;
		numExecutions = 1;
	} else {
		startVertex = rand() % numVertices;
	}
    // Executar Dijkstra para cada vértice
    for (int i = 0; i < numExecutions; ++i) {
		// Imprimir os caminhos mínimos em arquivo, se a flag de depuração estiver ativada.
        if (debugFlag) {
			pair<vector<int>, vector<int>> result = graph.dijkstra(startVertex);
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
				for (int targetVertex = 2; targetVertex <= numVertices; targetVertex++) {
					vector<int> path = graph.reconstructPath(startVertex, targetVertex - 1, result.second);
					// Exibe o caminho mínimo no formato especificado
					spath << "[" << startVertex + 1 << "," << targetVertex << "](" << result.first[targetVertex - 1] << ") ";
					for (int i = 0; i < path.size(); ++i) {
						spath << path[i] << " ";
					}
					spath << endl;
				}
			}
			spath.close();
        } else {		
			cout << "\nstartVertex = " << (startVertex + 1) << endl;
			pair<vector<int>, vector<int>> result = graph.dijkstra(startVertex);

			for (int targetVertex = 1; targetVertex <= numVertices; targetVertex++) {
				vector<int> path = graph.reconstructPath(startVertex, targetVertex - 1, result.second);
				// Exibe o caminho mínimo no formato especificado
				cout << "\n[" << startVertex + 1 << "," << targetVertex << "](" << result.first[targetVertex - 1] << ") ";
				for (int i = 0; i < path.size(); ++i) {
					cout << path[i] << " ";
				}
			}				        
			cout << endl;
			startVertex = rand() % numVertices;
		}
    }
    return 0;
}

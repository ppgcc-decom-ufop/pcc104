/**
 * PCC104 - Projeto e Analise de Algoritmos
 * Departamento de Computacao - Universidade Federal de Ouro Preto - MG
 * Professor: Pedro Silva
 * Aluno: Fernando dos Santos Alves Fernandes
 * Trabalho Pratico III (TP III) - Tarefa 2: Problema do Minimo Numero de Cliques
 * Algoritmos aproximados: Hill-climbing, Simulated Annealing, Tabu Search
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
#include <set>
#include <numeric>
#include <unordered_set>
#include <unordered_map>

using namespace std;

class Graph {
public:
    int V;
    vector<vector<int> > adj;
	double maxTime;
    chrono::time_point<chrono::_V2::system_clock, chrono::duration<long long int, ratio<1ll, 1000000000ll> > > initTime;
    chrono::time_point<chrono::_V2::system_clock, chrono::duration<long long int, ratio<1ll, 1000000000ll> > > endTime;
    chrono::duration<long long int, ratio<1ll, 1000000000ll> > elapsedTime;	

    Graph(int V) : V(V), adj(V) {}

    void addEdge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    bool isClique(const vector<int>& vertices) {
        for (size_t i = 0; i < vertices.size(); ++i) {
            for (size_t j = i + 1; j < vertices.size(); ++j) {
                if (find(adj[vertices[i]].begin(), adj[vertices[i]].end(), vertices[j]) == adj[vertices[i]].end()) {
                    return false;
                }
            }
        }
        return true;
    }

    vector<int> getRandomClique(const unordered_set<int>& coveredVertices) {
        vector<int> clique;
        int start;
        do {
            start = rand() % V;
        } while (coveredVertices.find(start) != coveredVertices.end());
        clique.push_back(start);
        return clique;
    }

    vector<int> getNeighbors(const vector<int>& clique, const unordered_set<int>& coveredVertices) {
        vector<int> neighbors;
        for (int v : clique) {
            for (int u : adj[v]) {
                if (find(clique.begin(), clique.end(), u) == clique.end() && coveredVertices.find(u) == coveredVertices.end()) {
                    neighbors.push_back(u);
                }
            }
        }
        return neighbors;
    }
	
	vector<int> hillClimbing(const unordered_set<int>& coveredVertices) {
		// Inicializa a clique atual com um vertice aleatorio nao coberto
		vector<int> current = getRandomClique(coveredVertices);

		// Loop principal
		while (true) {
			// Obtem os vizinhos da clique atual que nao estao cobertos
			vector<int> neighbors = getNeighbors(current, coveredVertices);

			// Inicializa a melhor vizinhanca como a clique atual
			vector<int> bestNeighbor = current;

			// Itera sobre todos os vizinhos possiveis
			for (int neighbor : neighbors) {
				// Cria uma nova clique adicionando o vizinho ao atual
				vector<int> newClique = current;
				newClique.push_back(neighbor);

				// Verifica se a nova clique e de fato uma clique e se e maior que a melhor vizinhanca atual
				if (isClique(newClique) && newClique.size() > bestNeighbor.size()) {
					bestNeighbor = newClique; // Atualiza a melhor vizinhanca encontrada
				}
			}

			// Se nao houve melhoria na clique atual, termina o algoritmo
			if (bestNeighbor.size() == current.size()) {
				break;
			}

			// Atualiza a clique atual para a melhor vizinhanca encontrada
			current = bestNeighbor;
		}

		// Retorna a melhor clique encontrada
		return current;
	}
	
	vector<int> simulatedAnnealing(const unordered_set<int>& coveredVertices, int maxIterations, double initialTemp, double coolingRate) {
		// Inicializa a clique atual e a melhor clique com um vertice aleatorio nao coberto
		vector<int> current = getRandomClique(coveredVertices);
		vector<int> best = current;

		// Inicializa a temperatura
		double temp = initialTemp;

		// Loop principal
		for (int iteration = 0; iteration < maxIterations; ++iteration) {
			// Obtem os vizinhos da clique atual que nao estao cobertos
			vector<int> neighbors = getNeighbors(current, coveredVertices);

			// Se nao ha vizinhos disponiveis, encerra o algoritmo
			if (neighbors.empty()) break;

			// Seleciona aleatoriamente um vizinho para explorar
			int neighbor = neighbors[rand() % neighbors.size()];
			vector<int> newClique = current;
			newClique.push_back(neighbor);

			// Verifica se a nova clique e de fato uma clique
			if (isClique(newClique)) {
				// Calcula a variacao de energia (tamanho da nova clique - tamanho da clique atual)
				double deltaE = newClique.size() - current.size();

				// Aceita a nova clique de acordo com uma probabilidade baseada na temperatura atual
				if (deltaE > 0 || exp(deltaE / temp) > ((double) rand() / RAND_MAX)) {
					current = newClique; // Aceita a nova clique
				}

				// Atualiza a melhor clique se a atual for melhor
				if (current.size() > best.size()) {
					best = current;
				}
			}

			// Reduz a temperatura de acordo com a taxa de resfriamento
			temp *= coolingRate;
		}

		// Retorna a melhor clique encontrada
		return best;
	}
	
	vector<int> tabuSearch(const unordered_set<int>& coveredVertices, int maxIterations, int tabuTenure) {
		// Inicializa a clique atual e a melhor clique com um vertice aleatorio
		vector<int> current = getRandomClique(coveredVertices);
		vector<int> best = current;

		// Lista tabu para armazenar os cliques visitados recentemente
		deque<vector<int>> tabuList;
		tabuList.push_back(current); // Adiciona a clique atual na lista tabu

		// Loop principal da Busca Tabu
		for (int iteration = 0; iteration < maxIterations; ++iteration) {
			// Obtem os vizinhos do clique atual que nao estao cobertos
			vector<int> neighbors = getNeighbors(current, coveredVertices);

			// Se nao ha vizinhos disponiveis, encerra o algoritmo
			if (neighbors.empty()) break;

			// Inicializa a melhor vizinhanca com a clique atual
			vector<int> bestNeighbor = current;

			// Itera sobre todos os vizinhos possiveis
			for (int neighbor : neighbors) {
				// Cria uma nova clique adicionando o vizinho na clique atual
				vector<int> newClique = current;
				newClique.push_back(neighbor);

				// Verifica se a nova clique e um clique e se nao esta na lista tabu
				if (isClique(newClique) && (find(tabuList.begin(), tabuList.end(), newClique) == tabuList.end())) {
					bestNeighbor = newClique; // Atualiza a melhor vizinhanca encontrada
					break; // Interrompe a busca apos encontrar uma vizinhanca valida
				}
			}

			// Se a melhor vizinhanca encontrada for melhor que a clique atual
			if (bestNeighbor.size() > current.size()) {
				current = bestNeighbor; // Atualiza o clique atual para a melhor vizinhanca encontrada

				// Atualiza a melhor clique global, se o clique atual for melhor
				if (current.size() > best.size()) {
					best = current;
				}

				// Adiciona o clique atual na lista tabu
				tabuList.push_back(current);

				// Remove o clique mais antigo da lista tabu se a lista atingir o tamanho maximo permitido
				if (tabuList.size() > tabuTenure) {
					tabuList.pop_front();
				}
			} else {
				break; // Encerra o loop se nao houver melhoria no clique atual
			}
		}

		// Retorna o melhor clique encontrado
		return best;
	}
	
    bool isSubset(const vector<int>& clique, const vector<vector<int>>& cliques) {
        for (const auto& existingClique : cliques) {
            unordered_set<int> existingSet(existingClique.begin(), existingClique.end());
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

	vector<vector<int>> vertexCoverByCliques(int timeLimit, string strategy) {
		// Inicializa o contador de tempo e a lista de cliques encontrados
		auto start = chrono::steady_clock::now();
		vector<vector<int>> cliques;
		unordered_set<int> coveredVertices;
		int maxIterations = 1000; // Numero maximo de iteracoes para as estrategias Simulated Annealing e Tabu Search
		double initialTemp = 1000.0; // Temperatura inicial para o Simulated Annealing
		double coolingRate = 0.95; // Taxa de resfriamento para o Simulated Annealing
		int tabuTenure = 10; // Tamanho da lista tabu para o Tabu Search

		// Loop principal para encontrar cobertura minima por cliques dentro do limite de tempo
		while (coveredVertices.size() < V) {
			vector<int> clique;

			// Seleciona a estrategia de busca de acordo com o parâmetro 'strategy'
			if (strategy == "hc") // Hill-climbing
				clique = hillClimbing(coveredVertices);
			else if (strategy == "sa") // Simulated Annealing
				clique = simulatedAnnealing(coveredVertices, maxIterations, initialTemp, coolingRate);
			else if (strategy == "ts") // Tabu Search
				clique = tabuSearch(coveredVertices, maxIterations, tabuTenure);

			// Verifica se o clique encontrado nao e um subconjunto de cliques ja encontrados
			if (!isSubset(clique, cliques)) {
				cliques.push_back(clique); // Adiciona o clique aa lista de cliques encontrados

				// Marca os vertices cobertos pela nova clique na cobertura
				for (int v : clique) {
					coveredVertices.insert(v);
				}
			}

			// Verifica o tempo decorrido e encerra o loop se o limite de tempo definido for atingido
			auto end = chrono::steady_clock::now();
			auto duration = chrono::duration_cast<chrono::seconds>(end - start).count();
			if (duration >= timeLimit) {
				break;
			}
		}

		// Ordena os cliques encontrados do maior para o menor tamanho
		sort(cliques.begin(), cliques.end(), [](const vector<int>& a, const vector<int>& b) {
			return a.size() > b.size();
		});

		// Retorna a lista ordenada de cliques encontrados
		return cliques;
	}
};

int main(int argc, char* argv[]) {
    if (argc != 4) {
        cerr << "Uso: " << argv[0] << " <arquivoProblema> <tempoMaximo> <strategy>" << endl;
        return 1;
    }

    string filename = argv[1];
    int timetorun = atoi(argv[2]);
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
	
	vector<vector<int>> cliques;
	cliques = graph.vertexCoverByCliques(timetorun, strategy);

    cout << cliques.size() << endl;
    for (const auto& clique : cliques) {
        for (int v : clique) {
            cout << v + 1 << " ";
        }
        cout << endl;
    }

    return 0;
}

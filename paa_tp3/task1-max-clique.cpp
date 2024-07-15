/**
 * PCC104 - Projeto e Analise de Algoritmos
 * Departamento de Computacao - Universidade Federal de Ouro Preto - MG
 * Professor: Pedro Silva
 * Aluno: Fernando dos Santos Alves Fernandes
 * Trabalho Pratico III (TP III) - Tarefa 1: Problema do Clique Maximo
 * Algoritmos: Hill Climbing, Simulated Annealing e Tabu Search
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
#include <cmath>
#include <set>
#include <numeric>
#include <unordered_set>

using namespace std;

class Graph {
private:
    int V;
    vector<vector<int> > adj;
	double maxTime;
    chrono::time_point<chrono::_V2::system_clock, chrono::duration<long long int, ratio<1ll, 1000000000ll> > > initTime;
    chrono::time_point<chrono::_V2::system_clock, chrono::duration<long long int, ratio<1ll, 1000000000ll> > > endTime;
    chrono::duration<long long int, ratio<1ll, 1000000000ll> > elapsedTime;
	
public:
    Graph(int V) : V(V) {
        adj.resize(V);
    }
	
    void setInitialAndMaxTime(chrono::time_point<chrono::_V2::system_clock, chrono::duration<long long int, ratio<1ll, 1000000000ll> > > initial, double time) {
        initTime = initial;
        maxTime = time;
    }

    auto getElasedTime(){
        return elapsedTime.count();
    }
	
	void printElapsedTime(){
		cout << "Tempo de execucao: " << fixed << setprecision(4) << double (elapsedTime.count() * 1e-9) << " s."<< endl;
	}	

    void addEdge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    bool isClique(const vector<int>& vertices) { // Avaliacao se os vertices formam um clique
        for (size_t i = 0; i < vertices.size(); ++i) {
            for (size_t j = i + 1; j < vertices.size(); ++j) {
                if (find(adj[vertices[i]].begin(), adj[vertices[i]].end(), vertices[j]) == adj[vertices[i]].end()) {
                    return false;
                }
            }
        }
        return true;
    }

    vector<int> getRandomClique() { // Estado inicial
        vector<int> clique;
        int start = rand() % V;
        clique.push_back(start);
        return clique;
    }

    vector<int> getNeighbors(const vector<int>& clique) { // Busca pela vizinhanca
        vector<int> neighbors;
        for (int v : clique) {
            for (int u : adj[v]) {
                if (find(clique.begin(), clique.end(), u) == clique.end()) {
                    neighbors.push_back(u);
                }
            }
        }
        return neighbors;
    }

    vector<int> hillClimbing() {
		// Verifica quanto tempo decorrido
		endTime = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::nanoseconds>(endTime - initTime);
		
		// Inicializa uma clique aleatoria
        vector<int> current = getRandomClique();
		
		// Loop principal que continua ate que nao se possa encontrar um clique maior ou o tempo maximo seja atingido
        while (true && (elapsed < chrono::duration_cast<chrono::nanoseconds>(chrono::duration<double>(maxTime)))) {
            // Obtem os vizinhos da clique atual
			vector<int> neighbors = getNeighbors(current);
			// Define a melhor clique vizinha como a clique atual
            vector<int> bestNeighbor = current;
			
			// Itera sobre os vizinhos
            for (int neighbor : neighbors) {
				// Cria um novo clique adicionando o vizinho na clique atual
                vector<int> newClique = current;
                newClique.push_back(neighbor);
				
				// Verifica se o novo conjunto e uma clique e se e maior que a melhor vizinha atual
                if (isClique(newClique) && newClique.size() > bestNeighbor.size()) {
                    bestNeighbor = newClique;
                }
            }
			
			// Se a melhor vizinha tem o mesmo tamanho da clique atual, interrompe o loop
            if (bestNeighbor.size() == current.size()) {
                break;
            }
			// Atualiza a clique atual para a melhor vizinha encontrada
            current = bestNeighbor;
			
			// Atualiza o tempo decorrido
			endTime = chrono::high_resolution_clock::now();
			auto elapsed = chrono::duration_cast<chrono::nanoseconds>(endTime - initTime);
			elapsedTime = elapsed;
        }
		// Retorna a clique atual, que e a maior clique encontrada
        return current;
    }
	
	// Funcao para calcular a probabilidade de aceitacao de uma solucao pior
    double probability(double delta, double temperature) {
        return exp(delta / temperature);
    }

    vector<int> simulatedAnnealing(double initialTemperature, double coolingRate) {
		// Verifica quanto tempo decorrido
		endTime = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::nanoseconds>(endTime - initTime);
		
        // Inicializa uma clique aleatoria e define o melhor clique inicial como esta clique
		vector<int> current = getRandomClique();
        vector<int> best = current;
		// Define a temperatura inicial
        double temperature = initialTemperature;

        // Loop principal que continua ate a temperatura cair abaixo de 1 ou o tempo maximo ser atingido
		while (temperature > 1 && (elapsed < chrono::duration_cast<chrono::nanoseconds>(chrono::duration<double>(maxTime)))) {
			// Obtem os vizinhos da clique atual
            vector<int> neighbors = getNeighbors(current);
			// Se nao houver vizinhos, interrompe o loop
            if (neighbors.empty()) {
                break;
            }
			// Seleciona um vizinho aleatorio
            int idx = rand() % neighbors.size();
            vector<int> newClique = current;
            newClique.push_back(neighbors[idx]);

            // Verifica se o novo conjunto e uma clique
			if (isClique(newClique)) {
				// Aceita a nova clique se for maior que a clique atual ou com uma certa probabilidade mesmo que seja for menor
                if (newClique.size() > current.size() || probability(newClique.size() - current.size(), temperature) > (double)rand() / RAND_MAX) {
                    current = newClique;
                }
				// Atualiza o melhor clique encontrado ate agora
                if (current.size() > best.size()) {
                    best = current;
                }
            }
			// Reduz a temperatura
            temperature *= coolingRate;
			
			// Atualiza o tempo decorrido
			endTime = chrono::high_resolution_clock::now();
			auto elapsed = chrono::duration_cast<chrono::nanoseconds>(endTime - initTime);
			elapsedTime = elapsed;
        }
		// Retorna a melhor clique encontrada
        return best;
    }
	
    vector<int> tabuSearch(int maxIterations, int maxTabuSize) {
		// Verifica quanto tempo decorrido
		endTime = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::nanoseconds>(endTime - initTime);
		
		// Inicializa uma clique aleatoria e define o melhor clique inicial como esta clique
        vector<int> current = getRandomClique();
        vector<int> best = current;
		
		// Define uma lista tabu para evitar ciclos
        unordered_set<int> tabuList;
        int iterations = 0;
		
		// Loop principal que continua ate o maximo de iteracoes ou tempo maximo
        while (iterations < maxIterations && (elapsed < chrono::duration_cast<chrono::nanoseconds>(chrono::duration<double>(maxTime)))) {
			// Obtem os vizinhos da clique atual
            vector<int> neighbors = getNeighbors(current);
			// Se nao houver vizinhos, interrompe o loop
            if (neighbors.empty()) {
                break;
            }

            // Variavel para armazenar o melhor vizinho encontrado
			vector<int> bestNeighbor;
			// Itera sobre os vizinhos
            for (int neighbor : neighbors) {
				// Se o vizinho nao esta na lista tabu
                if (tabuList.find(neighbor) == tabuList.end()) {
					// Cria um novo clique adicionando o vizinho na clique atual
                    vector<int> newClique = current;
                    newClique.push_back(neighbor);

                    // Verifica se o novo conjunto e uma clique e se e maior que o melhor vizinho atual
					if (isClique(newClique) && newClique.size() > bestNeighbor.size()) {
                        bestNeighbor = newClique;
                    }
                }
            }

            // Se nenhum vizinho valido foi encontrado, interrompe o loop
			if (bestNeighbor.empty()) {
                break;
            }

            // Atualiza a clique atual para o melhor vizinho encontrado
			current = bestNeighbor;

            // Se a clique atual e maior que a melhor encontrada ate agora, atualiza a melhor clique e limpa a lista tabu
			if (current.size() > best.size()) {
                best = current;
                tabuList.clear();
            }

            // Adiciona o ultimo vertice da clique atual na lista tabu
			tabuList.insert(current.back());
			// Se a lista tabu excede o tamanho maximo permitido, remove o elemento mais antigo
            if (tabuList.size() > maxTabuSize) {
                tabuList.erase(tabuList.begin());
            }

            // Incrementa o contador de iteracoes e atualiza o tempo decorrido
			iterations++;
			
			// Contabiliza o tempo de execucao ate o momento
			endTime = chrono::high_resolution_clock::now();
			auto elapsed = chrono::duration_cast<chrono::nanoseconds>(endTime - initTime);
			elapsedTime = elapsed;
        }
		// Retorna a melhor clique encontrada
        return best;
    }	
};

int main(int argc, char* argv[]) {
    if (argc != 4) {
        cerr << "Uso: " << argv[0] << " <arquivoProblema> <tempoMaximo> <strategy>" << endl;
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
	
	srand(time(0));
	
	auto initTime = chrono::high_resolution_clock::now();
    graph.setInitialAndMaxTime(initTime, timetorun);
	
	vector<int> clique;
	if (strategy == "hc") // Hill-climbing
		clique = graph.hillClimbing();
	else if (strategy == "sa") // Simulated Annealing
		clique = graph.simulatedAnnealing(1000, 0.95); // Temperatura inicial e temperatura minima
	else if (strategy == "ts") // Tabu Search
		clique = graph.tabuSearch(10000, 7); // Maximo de iteracoes e tamanho maximo da lista tabu
	graph.printElapsedTime();

	cout << clique.size() << endl;
    for (int v : clique) {
        cout << v + 1 << " ";
    }
    cout << endl;
    return 0;
}

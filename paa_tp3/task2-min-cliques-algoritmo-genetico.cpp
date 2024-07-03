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
#include <random>
#include <set>
#include <map>
#include <utility>
#include <ctime>
#include <cstdlib> // srand, rand

using namespace std;
using namespace std::chrono;

class Graph {
public:
    int V;
    vector<vector<int>> adj;

    Graph(int V) : V(V), adj(V) {}

    void addEdge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    vector<int> neighbors(int node) const {
        return adj[node];
    }
};

using Individual = vector<int>;

Individual generate_individual(int V) {
    Individual individual(V);
    iota(individual.begin(), individual.end(), 0);
    shuffle(individual.begin(), individual.end(), default_random_engine(rand()));
    return individual;
}

vector<set<int>> evaluate_individual(const Individual &individual, const Graph &graph) {
    set<int> covered_vertices;
    vector<set<int>> cliques;

    for (int node : individual) {
        if (covered_vertices.find(node) == covered_vertices.end()) {
            set<int> clique = {node};
            for (int neighbor : graph.neighbors(node)) {
                if (covered_vertices.find(neighbor) == covered_vertices.end()) {
                    set<int> potential_clique = clique;
                    potential_clique.insert(neighbor);
                    bool is_clique = true;
                    for (int a : potential_clique) {
                        for (int b : potential_clique) {
                            if (a != b && find(graph.adj[a].begin(), graph.adj[a].end(), b) == graph.adj[a].end()) {
                                is_clique = false;
                                break;
                            }
                        }
                        if (!is_clique) break;
                    }
                    if (is_clique) {
                        clique = potential_clique;
                    }
                }
            }
            cliques.push_back(clique);
            covered_vertices.insert(clique.begin(), clique.end());
            if (covered_vertices.size() == graph.V) {
                break;
            }
        }
    }

    if (covered_vertices.size() == graph.V) {
        return cliques;
    } else {
        return {};
    }
}

void mutate_individual(Individual &individual) {
    int idx1 = rand() % individual.size();
    int idx2 = rand() % individual.size();
    swap(individual[idx1], individual[idx2]);
}

void crossover(Individual &ind1, Individual &ind2) {
    int size = ind1.size();
    int cxpoint1 = rand() % size;
    int cxpoint2 = rand() % size;
    if (cxpoint2 >= cxpoint1) {
        cxpoint2++;
    } else {
        swap(cxpoint1, cxpoint2);
    }
    for (int i = cxpoint1; i < cxpoint2; i++) {
        swap(ind1[i], ind2[i]);
    }
}

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
    srand(time(0));

    // Parameters
    int population_size = 50;
    int num_generations = 50;
    double crossover_prob = 0.7;
    double mutation_prob = 0.2;

    // Initialize population
    vector<Individual> population(population_size);
    generate(population.begin(), population.end(), [&]() { return generate_individual(graph.V); });

    auto start = high_resolution_clock::now();

    for (int gen = 0; gen < num_generations; gen++) {
        // Check time limit
        auto now = high_resolution_clock::now();
        duration<double> elapsed = now - start;
        if (elapsed.count() >= timetorun) {
            break;
        }

        // Evaluate individuals
        vector<pair<int, Individual>> fitness_population;
        for (auto &ind : population) {
            auto cliques = evaluate_individual(ind, graph);
            int fitness = cliques.size();
            fitness_population.push_back({fitness, ind});
        }

        // Sort by fitness
        sort(fitness_population.begin(), fitness_population.end());

        // Selection
        vector<Individual> new_population;
        for (int i = 0; i < population_size; i++) {
            new_population.push_back(fitness_population[i].second);
        }

        // Crossover
        for (int i = 0; i < population_size; i += 2) {
            if ((double)rand() / RAND_MAX < crossover_prob && i + 1 < population_size) {
                crossover(new_population[i], new_population[i + 1]);
            }
        }

        // Mutation
        for (int i = 0; i < population_size; i++) {
            if ((double)rand() / RAND_MAX < mutation_prob) {
                mutate_individual(new_population[i]);
            }
        }

        // Replace population
        population = new_population;

        // Print statistics
        int min_fit = fitness_population.front().first;
        int max_fit = fitness_population.back().first;
        double mean_fit = 0;
        for (auto &fp : fitness_population) {
            mean_fit += fp.first;
        }
        mean_fit /= population_size;

        // Uncomment the following line to print generation statistics
        // cout << "Generation " << gen << ": Min " << min_fit << ", Max " << max_fit << ", Avg " << mean_fit << endl;
    }

    // Find the best individual
    auto best_individual = *max_element(population.begin(), population.end(), [&](const Individual &a, const Individual &b) {
        return evaluate_individual(a, graph).size() < evaluate_individual(b, graph).size();
    });
    auto best_cliques = evaluate_individual(best_individual, graph);

    cout << best_cliques.size() << endl;
    for (auto &clique : best_cliques) {
        cout << "";
        for (int node : clique) {
            cout << node + 1 << " "; // Adjusted to print 1-based node indices
        }
        cout << "\n";
    }
    cout << endl;

    return 0;
}

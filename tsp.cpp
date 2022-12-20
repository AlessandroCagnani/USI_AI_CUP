#include <iostream>
#include <fstream>
#include <math.h>
#include <ctime>
#include <string>
#include <vector>
#include <iomanip>
#include <random>
#include <algorithm>
#include <climits>
#include <cmath>
#include <numeric>

#include "point.hpp"

using namespace std;

/**
 * @brief TSP problem class
 *
 * @param points list of euclidean points for this problem
 *
 */
class problem
{
public:
    string name;
    string type;
    string comment;
    int dimension;
    string edge_weight_type;
    int best_known;
    vector<point> points;
    vector<vector<int>> distance_matrix;

    problem(const char *problem)
    {
        string path = "problems/";
        path.append(problem);

        ifstream problem_file;
        problem_file.open(path.c_str());

        if (!problem_file.is_open())
        {
            cout << "problem file not found" << endl;
            abort();
        }
        string line;
        // state 0: fetch info, state 1: fetch points
        int state = 0;
        int line_num = 0;
        while (getline(problem_file, line))
        {
            if (line == "NODE_COORD_SECTION")
            {
                state = 1;
                continue;
            }

            if (state == 0)
            {
                if (line.find("NAME") != string::npos)
                {
                    name = line.substr(line.find(":") + 1);
                }
                else if (line.find("TYPE") != string::npos && line_num == 1)
                {
                    type = line.substr(line.find(":") + 1);
                }
                else if (line.find("COMMENT") != string::npos)
                {
                    comment = line.substr(line.find(":") + 1);
                }
                else if (line.find("DIMENSION") != string::npos)
                {
                    dimension = stoi(line.substr(line.find(":") + 1));
                }
                else if (line.find("EDGE_WEIGHT_TYPE") != string::npos)
                {
                    edge_weight_type = line.substr(line.find(":") + 1);
                }
                else if (line.find("BEST_KNOWN") != string::npos)
                {
                    best_known = stoi(line.substr(line.find(":") + 1));
                }
            }
            else if (state == 1)
            {

                if (line == "EOF")
                {
                    break;
                }
                point p;

                int space_index = line.find(" ");
                int second_space_index = line.substr(space_index + 1).find(" ");
                // cout << setprecision(15) << "cord x:" << stod(line.substr(space_index + 1, second_space_index)) << "{X}" << endl;
                // cout << setprecision(15) << "cord y:" << stod(line.substr(space_index + 1 + second_space_index + 1)) << "{Y}" << endl;
                // int index;
                // problem_file >> index >> p.coord_x >> p.coord_y;
                // cout << "index: " << index << " x: " << p.coord_x << " y: " << p.coord_y << endl;
                p.coord_x = stod(line.substr(space_index + 1, second_space_index));
                p.coord_y = stod(line.substr(space_index + 1 + second_space_index + 1));
                points.push_back(p);
            }
            line_num++;
        }
        problem_file.close();

        distance_matrix = vector<vector<int>>(dimension, vector<int>(dimension));
        for (int i = 0; i < dimension; i++)
        {
            for (int j = 0; j < dimension; j++)
            {
                if (i == j)
                {
                    double distance = INFINITY;
                    distance_matrix[i][j] = distance;
                }
                else
                {

                    int distance = euclidean_distance(points[i], points[j]);
                    if (distance == 0)
                    {
                        cout << "distance between " << i << " and " << j << " is " << distance << endl;
                        cout << "point " << i << " x: " << points[i].coord_x << " y: " << points[i].coord_y << endl;
                        cout << "point " << j << " x: " << points[j].coord_x << " y: " << points[j].coord_y << endl;
                    }
                    distance_matrix[i][j] = distance;
                }
            }
        }
        // print distance matrix
        // for (int i = 0; i < dimension; i++)
        // {
        //     for (int j = 0; j < dimension; j++)
        //     {
        //         cout << distance_matrix[i][j] << " ";
        //     }
        //     cout << endl;
        // }
    }

    vector<int> nearest_neighbour(int starting_node = 0)
    {
        vector<int> path(dimension);
        vector<char> visited(dimension, 0);

        path[0] = starting_node;
        visited[0] = 1;

        for (int c = 1; c < dimension; c++)
        {
            int current = path[c - 1];
            int next_city;
            int dist = INT_MAX;
            for (int i = 0; i < dimension; i++)
            {
                if (visited[i] == 0 && distance_matrix[current][i] < dist)
                {
                    next_city = i;
                    dist = distance_matrix[current][i];
                }
            }
            path[c] = next_city;
            visited[next_city] = 1;
        }

        return path;
    }

    vector<int> nearest_neighbors_algorithm(vector<bool> &already_visited_nodes, int starting_node)
    {
        vector<int> solution;
        solution.push_back(starting_node);
        already_visited_nodes[starting_node] = true;
        int current_node = starting_node;
        for (int i = 0; i < dimension - 1; i++)
        {
            int nearest_neighbor = -1;
            float min_distance = INFINITY;
            for (int j = 0; j < dimension; j++)
            {
                if (already_visited_nodes[j] == false)
                {
                    if (distance_matrix[current_node][j] < min_distance)
                    {
                        nearest_neighbor = j;
                        min_distance = distance_matrix[current_node][j];
                    }
                }
            }
            solution.push_back(nearest_neighbor);
            already_visited_nodes[nearest_neighbor] = true;
            current_node = nearest_neighbor;
        }
        return solution;
    }

    vector<int> best_nearest_neighbors()
    {
        vector<int> solution;
        int n = dimension;
        float min_length = INFINITY;
        for (int i = 0; i < n; i++)
        {
            vector<bool> already_visited_nodes(n, false);
            // solution returned by the nearest neighbors algorithm
            vector<int> nn = nearest_neighbors_algorithm(already_visited_nodes, i);
            int current_length = get_cost(nn);
            if (current_length < min_length)
            {
                solution = nn;
                min_length = current_length;
            }
        }
        return solution;
    }

    double get_cost(vector<int> solution)
    {
        double cost = 0;
        for (int i = 0; i < dimension - 1; i++)
        {
            cost += distance_matrix[solution[i]][solution[i + 1]];
        }
        cost += distance_matrix[solution[dimension - 1]][solution[0]];
        return cost;
    }

    bool valid(const vector<int> &solution) const
    {
        vector<char> visited(dimension, 0);
        for (int i : solution)
        {
            if (visited[i])
                return false;
            visited[i] = 1;
        }
        return true;
    }

    void print_fields()
    {
        cout << "name: " << name << endl;
        cout << "type: " << type << endl;
        cout << "comment: " << comment << endl;
        cout << "dimension: " << dimension << endl;
        cout << "edge_weight_type: " << edge_weight_type << endl;
        cout << "best_known: " << best_known << endl;
        cout << "points: " << endl;
        for (int i = 0; i < points.size(); i++)
        {
            cout << points[i].coord_x << " " << points[i].coord_y << endl;
        }
        cout << "distance_matrix: " << endl;
        for (int i = 0; i < dimension; i++)
        {
            cout << "[";
            for (int j = 0; j < dimension; j++)
            {
                cout << distance_matrix[i][j] << " ";
            }
            cout << "]" << endl;
        }
    }

    vector<int> random_solution()
    {
        vector<int> solution;
        for (int i = 0; i < dimension; i++)
        {
            solution.push_back(i);
        }
        random_shuffle(solution.begin(), solution.end());
        return solution;
    }
};

// random neighbors with 2opt
vector<int> random_neighbors(vector<int> solution)
{
    int size = solution.size();
    int random_index1 = rand() % size;
    int random_index2 = rand() % size;
    while (random_index1 == random_index2)
    {
        random_index2 = rand() % size;
    }
    if (random_index1 > random_index2)
    {
        int temp = random_index1;
        random_index1 = random_index2;
        random_index2 = temp;
    }
    vector<int> new_solution = solution;
    for (int i = 0; i < (random_index2 - random_index1 + 1) / 2; i++)
    {
        int temp = new_solution[random_index1 + i];
        new_solution[random_index1 + i] = new_solution[random_index2 - i];
        new_solution[random_index2 - i] = temp;
    }
    return new_solution;
}

pair<vector<int>, bool> step2opt(vector<int> tsp_sequence, vector<vector<int>> matrix_dist)
{
    vector<int> new_tsp_sequence = tsp_sequence;
    bool reversed = false;
    int minimum_gain = 0;
    int best_i = -1;
    int best_j = -1;
    for (int i = 1; i < new_tsp_sequence.size() - 2; i++)
    {
        for (int j = i + 1; j < new_tsp_sequence.size() - 1; j++)
        {
            int old_links_len = matrix_dist[new_tsp_sequence[i]][new_tsp_sequence[i - 1]] + matrix_dist[new_tsp_sequence[j]][new_tsp_sequence[j + 1]];
            int changed_links_len = matrix_dist[new_tsp_sequence[j]][new_tsp_sequence[i - 1]] + matrix_dist[new_tsp_sequence[i]][new_tsp_sequence[j + 1]];
            // having a function to compute the gain value slowed down the computation of the solution
            int gain = changed_links_len - old_links_len;
            if (gain < minimum_gain)
            {
                reversed = true;
                best_i = i;
                best_j = j;
                minimum_gain = gain;
            }
        }
    }
    if (reversed)
    {
        reverse(new_tsp_sequence.begin() + best_i, new_tsp_sequence.begin() + best_j + 1);
    }
    return make_pair(new_tsp_sequence, reversed);
}

vector<int> loop2opt(vector<int> solution, problem *instance)
{
    bool reversed = true;
    vector<int> new_solution = solution;
    while (reversed)
    {
        pair<vector<int>, bool> step_solution = step2opt(new_solution, instance->distance_matrix);
        new_solution = step_solution.first;
        reversed = step_solution.second;
    }
    return new_solution;
}

int cost_change(vector<vector<int>> *distance_matrix, int n1, int n2, int n3, int n4)
{
    return (*distance_matrix)[n1][n3] + (*distance_matrix)[n2][n4] - (*distance_matrix)[n1][n2] - (*distance_matrix)[n3][n4];
}


vector<int> random_2opt(std::mt19937 &gen, problem &p, vector<int> &current)
{
    std::uniform_real_distribution<> dis(0.0, 1.0);
    int n = current.size();

    int a = int(dis(gen) * n) % n;
    int b = (a + 1) % n;
    int d = int(dis(gen) * n) % n;
    int c = (d + 1) % n;

    vector<int> next = current;
    reverse(next.begin() + a, next.begin() + d + 1);

    return next;
}

void swap(int i, int j, vector<int> *solution)
{
    int temp = (*solution)[i];
    (*solution)[i] = (*solution)[j];
    (*solution)[j] = temp;
}



vector<int> annealing(problem p, int seed, double T_max = 100, double T_min = 1, double cooling = 0.99)
{
    // cout << "[annealing]" << endl;
    uniform_real_distribution<> dis(0.0, 1.0);
    mt19937 gen(seed);
    double T = T_max;

    // vector<int> solution = p.random_solution();
    vector<int> solution = p.best_nearest_neighbors();
    // vector<int> solution = p.nearest_neighbour();
    double current_cost = p.get_cost(solution);
    double best_cost = current_cost;
    vector<int> best_solution = solution;

    if (p.dimension > 440) {
        cooling = 0.9;
    }
    while (T > T_min)
    {
        // cout << "T: " << T << endl;
        int n = 10000;
        // cout << "dimension: " << p.dimension << endl;
        if (p.dimension > 440) {
            n = 10;
            // cout << "n: " << n << endl;
        }
        for (int i = 0; i < n; i++)
        {

            vector<int> new_solution;
            if (n == 10000) {
                new_solution = random_neighbors(solution);
            } else {
                // solution = loop2opt(solution, &p);
                new_solution = loop2opt(solution, &p);
                // cout << "ann->2opt" << endl;
            }
            // vector<int> new_solution1 = random_neighbors(solution);
            // vector<int> new_solution = loop2opt(solution, &p);
            // vector<int> new_solution1 = random_2opt(gen, p, solution);
            // vector<int> new_solution2 = two_dot_five_opt(solution, &p);

            // double new_cost1 = p.get_cost(new_solution1);
            // double new_cost2 = p.get_cost(new_solution2);
            
            // vector<int> new_solution;

            // if (new_cost1 < new_cost2) {
            //   new_solution = new_solution1;
            // } else {
            //   new_solution = new_solution2;
            // }

            double new_cost = p.get_cost(new_solution);
            double delta = new_cost - current_cost;

            if (delta < 0)
            {
                solution = new_solution;
                current_cost = new_cost;
                if (current_cost < best_cost)
                {
                    best_cost = current_cost;
                    best_solution = solution;
                    // cout << "[primo if] best cost: " << best_cost << endl;
                }
            }
            else
            {
                double prob = exp(-delta / T);
                if (dis(gen) < prob)
                {
                    solution = new_solution;
                    current_cost = new_cost;
                    // cout << "[secondo if] current cost: " << current_cost << endl;
                }
            }
        }
        T *= cooling;
    }
    return best_solution;
}


int main(int argc, const char *argv[])
{

    // take the time
    clock_t start, end;
    start = clock();

    long seed = rand();
    if (argc < 2)
    {
        cout << "pls run the program as follow ./tsp {tsp_problem name} {seed}" << endl;
        abort();
    }
    if (argc == 3)
    {
        seed = stoi(argv[2]);
    }

    problem p(argv[1]);

    vector<int> solution = annealing(p, seed);
    // vector<int> solution = p.best_nearest_neighbors();
    // solution = two_dot_five_opt(solution, &p);

    if (!p.valid(solution))
    {
        cout << "invalid solution" << endl;
    }
    else
    {
        cout << "---------------" << p.name << "--------------- \n"
             << endl;
        cout << "time: " << (double)(clock() - start) / CLOCKS_PER_SEC << "s" << endl;
        cout << "seed: " << seed << endl;

        cout << "\nbest solution: " << p.get_cost(solution) << endl;
        cout << "error: " << (double(p.get_cost(solution) - p.best_known) / p.best_known) * 100 << endl;
        cout << "solution: " << endl;
        for (int i = 0; i < solution.size(); i++)
        {
            cout << solution[i] << ", ";
        }
        cout << "\n"
             << endl;
    }

    return 0;
}

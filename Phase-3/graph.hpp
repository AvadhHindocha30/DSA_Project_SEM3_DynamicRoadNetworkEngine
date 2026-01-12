#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <vector>
#include <string>
#include <unordered_map>
#include <cmath>
#include <list>

struct Node {
    double lat, lon;
};

struct Edge {
    int id;
    int u, v;
    double length;
    double avg_time;
    bool oneWay;
    bool active;
    std::string road_type;
};

class Graph {
    std::vector<double> g_score;
    // std::vector<int> parent;
    std::vector<int> visited; 
    int current_token = 0;
    double find_heuristic(int u, int v);

public:
    std::vector<Node> nodes;
    std::unordered_map<int, Edge*> edges;
    std::vector<std::vector<Edge*>> adj;
    double max_speed;          // this needs to be used in the find_heuristic fn

    explicit Graph(int n);
    ~Graph();

    void addNode(int id, double lat, double lon);
    void addEdge(int id, int u, int v, double length, double avg_time, bool oneWay, const std::string &road_type);
    double euclidean_estimate(int u, int v);
    bool getShortestPath(int src, int dest, std::vector<int>& outPath, double& outTime);
};

#endif // GRAPH_HPP
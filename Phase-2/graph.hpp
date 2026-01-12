#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <limits>

static const double INF = std::numeric_limits<double>::infinity();

// Data structures for edges and nodes (same as Phase 1)
struct Edge {
    int id;                   // unique ID for every edge  
    int u,v;                  // vertices id between this edge is present
    double length;
    double average_time;

    bool oneWay;              // true if way is there from u to v only

    std::vector<double> speed_profile;
    std::string road_type;
    bool active; 
};

// node struct

struct Node {
    int id;
    double lat, lon;

    std::vector<std::string> pois;
};

//each path with its properties
struct Path {
    std::vector<int> nodes;
    std::vector<int> edges;
    double cost;
    double penalty;
};

class Graph {
public:

    int numNodes;
    std::vector<Node> nodes;
    std::unordered_map<int,Edge> edges;

    std::vector<std::vector<Edge*>> adj;

    Graph(int n);
    ~Graph();

    void addNode(int id, double lat, double lon, const std::vector<std::string>& pois);
    void addEdge(int id, int u, int v, double length, double average_time, bool oneway,
                 const std::string &road_type, const std::vector<double>& speed_profile);
    bool removeEdge(int edge_id);
    bool modifyEdge(int edge_id, const std::vector<double>& new_speed_profile,
                    double new_length, double new_avg_time, std::string newRoadType);

    // Shortest path (distance and time)
    bool shortestPathDistance(int src, int dest, std::vector<int>& outPath, double& outDist,
                                 const std::unordered_set<int>& forbidNodes,
                                 const std::unordered_set<int>& forbidEdges);
    bool shortestPathTime(int src, int dest, std::vector<int>& outPath, double& outTime,
                          const std::unordered_set<int>& forbidNodes = {},
                          const std::unordered_set<std::string>& forbidRoads = {});

    std::vector<Path> kShortestPathsExact(int src, int dest, int k, const std::string& mode);
    std::vector<Path> kShortestPathsHeuristic(int src, int dest, int k, double overlap_threshold);
    std::vector<std::tuple<int, int, double>> approxShortestDistances(const std::vector<std::pair<int,int>>& queries, double budget, double pct);
    
    double aStarDistance(int src, int dest, double pct);
};

#endif 

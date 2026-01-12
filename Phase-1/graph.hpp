#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>

// edge struct

struct Edge {
    int id;                   // unique ID for every edge  
    int u,v;                  // vertices id between this edge is present
    double length;
    double average_time;

    bool oneway;              // true if way is there from u to v only

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

// main graph class

class Graph {
public:

    int numNodes;
    std::vector<Node> nodes;
    std::unordered_map<int,Edge> edges;

    std::vector<std::vector<Edge*>> adj;


    Graph(int n);
    ~Graph();

    // graph_util functions

    void addNode(int id, double lat, double lon, const std::vector<std::string>& pois);
    void addEdge(int id, int u, int v, double length, double average_time, bool oneway,
                 const std::string &road_type, const std::vector<double>& speed_profile);
    bool removeEdge(int edge_id);
    bool modifyEdge(int edge_id, const std::vector<double>& new_speed_profile,
                    double new_length, double new_avg_time, std::string newRoadType);

   
    bool shortestPathDistance(int src, int dest, std::vector<int>& outPath, double& outDist,
                              const std::unordered_set<int>& forbidNodes = {},
                              const std::unordered_set<std::string>& forbidRoads = {});
    bool shortestPathTime(int src, int dest, std::vector<int>& outPath, double& outTime,
                          const std::unordered_set<int>& forbidNodes = {},
                          const std::unordered_set<std::string>& forbidRoads = {});

    std::vector<int> knnQueryEuclidean(const std::string& poiType, double queryLat,
                                       double queryLon, int k);
    std::vector<int> knnQueryShortestPath(const std::string& poiType, double queryLat,
                                          double queryLon, int k);


};
#endif 
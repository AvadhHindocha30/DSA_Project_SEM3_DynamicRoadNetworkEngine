#include "graph.hpp"
#include <algorithm>
#include <limits>
#include <queue>
#include <cmath>

#ifndef PI
#define PI 3.14159265358979323846
#endif
#define EPSILON 0.001         // Min possible time to cover an edge
#define EXPNASION_BUFFER 100
using pii = std::pair<double, int>; 


double deg_to_rad(double theta) {
    return theta * (PI / 180);
}

// Haversine Distance (Meters)
double haversine(const Node& a, const Node& b) {
    constexpr double R = 6371000.0;     // Earth's radius
    double dLat = deg_to_rad(b.lat - a.lat);
    double dLon = deg_to_rad(b.lon - a.lon);
    double lat1 = deg_to_rad(a.lat);
    double lat2 = deg_to_rad(b.lat);

    double x = std::sin(dLat / 2.0) * std::sin(dLat / 2.0) +
               std::cos(lat1) * std::cos(lat2) * std::sin(dLon / 2.0) * std::sin(dLon / 2.0);
    double theta = 2.0 * std::atan2(std::sqrt(x), std::sqrt(1.0 - x));
    return R * theta;      // This is the arc length on the Earth's which is assumed to be spherical
}


Graph::Graph(int n) {
    nodes.resize(n + 1);
    adj.resize(n + 1);
    g_score.resize(n + 1);
    visited.resize(n + 1);
    max_speed = 1.0;
}

Graph::~Graph() {
    for (auto& pair : edges) delete pair.second;
    // edges.clear();
}

void Graph::addNode(int id, double lat, double lon) {
    if (id >= nodes.size()) {
        nodes.resize(id + EXPNASION_BUFFER);
        adj.resize(id + EXPNASION_BUFFER);
        g_score.resize(id + EXPNASION_BUFFER);
        visited.resize(id + EXPNASION_BUFFER);
    }
    nodes[id] = {lat, lon};
}

void Graph::addEdge(int id, int u, int v, double length, double avg_time, bool oneWay, const std::string &road_type) {
    Edge* e = new Edge{id, u, v, length, avg_time, oneWay, true, road_type};
    edges[id] = e;
    adj[u].push_back(e);

    if (u >= adj.size() || v >= adj.size()) {
        int newSize = std::max(u, v) + EXPNASION_BUFFER;
        adj.resize(newSize);
        nodes.resize(newSize); 
        g_score.resize(newSize, std::numeric_limits<double>::max());
        visited.resize(newSize);
    }

    if (!oneWay) {
        Edge* rev = new Edge{id, v, u, length, avg_time, oneWay, true, road_type};
        adj[v].push_back(rev);
    }
    
    if (avg_time > EPSILON) {
        double speed = length / avg_time;
        max_speed = std::max(speed, max_speed);
    }
}

double Graph::find_heuristic(int u, int v) {
    return haversine(nodes[u], nodes[v]) / max_speed;
}

double Graph::euclidean_estimate(int u, int v) {
    if (u < 0 || u >= (int)nodes.size() || v < 0 || v >= (int)nodes.size()) return 1e9;
    return find_heuristic(u, v); 
}

bool Graph::getShortestPath(int src, int dest, std::vector<int>& outPath, double& outTime) {
    if (src == dest) {
        outPath = {src};
        outTime = 0.0;
        return true;
    }
    if (src >= g_score.size() || dest >= g_score.size()) return false;
    
    std::priority_queue<pii, std::vector<pii>, std::greater<pii>> pq;
    pq.push({find_heuristic(src, dest), src});
    g_score[src] = 0.0;
    visited[src] = ++current_token;
    std::vector<int> parent(nodes.size(), -1);
    bool found = false;

    // choose minimum value of f(n) where f(n) = g(n) + h(n) at node n
    // f(n) = total estimated trip time
    // g(n) = actual time taken to reach node 'n' from the start
    // h(n) = prediction by heuristic fn to get to the destination from node 'n'

    while (!pq.empty()) {
        double f = pq.top().first;
        int u = pq.top().second;
        pq.pop();
        if (visited[u] == current_token && f > g_score[u] + find_heuristic(u, dest)) continue;
        if (u == dest) { found = true; break; }

        for (Edge* e : adj[u]) {
            if (!e->active) continue;
            int v = e->v;
            double time = e->avg_time; 
            double g_v = (visited[v] == current_token) ? g_score[v] : std::numeric_limits<double>::max();
            if (g_score[u] + time < g_v) {
                g_score[v] = g_score[u] + time;
                parent[v] = u;
                visited[v] = current_token;
                pq.push({g_score[v] + find_heuristic(v, dest), v});
            }
        }
    }
    if (!found) return false; 

    outTime = g_score[dest];
    outPath.clear();
    int curr = dest;
    while (curr != -1) {
        outPath.push_back(curr);
        curr = parent[curr];
    }
    std::reverse(outPath.begin(), outPath.end());
    return true;
}
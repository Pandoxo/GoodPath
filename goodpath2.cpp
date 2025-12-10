#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <map>
#include <fstream>
#include <filesystem> 
#include <climits>
#include <random>
#include <stack>
#include <unordered_set>
#include "util.h"

namespace fs = std::filesystem;
using namespace std;

struct Edge {
    int to;
    int weight;
    
    Edge(int t, int w) : to(t), weight(w) {}
};

// Structure to store information about removed nodes
struct ReductionInfo {
    vector<vector<Edge>> graph;           // The reduced graph
    vector<bool> removed;                 // Which nodes were removed
    vector<pair<int, int>> replacedEdge;  // For each removed node: (neighbor1, neighbor2)
    map<pair<int, int>, vector<int>> edgeChains; // Map from edge to removed nodes
};

// Returns the reduced graph with removal information for path reconstruction
ReductionInfo reduceGraphWithInfo(const vector<vector<Edge>>& adjList) {
    int n = adjList.size();
    vector<vector<Edge>> graph = adjList;
    vector<bool> removed(n, false);
    vector<int> degree(n);
    vector<pair<int, int>> replacedEdge(n, {-1, -1});
    
    // Map from edge (u,v) to chain of removed nodes between them
    map<pair<int, int>, vector<int>> edgeToChain;
    
    // Initialize degrees and find initial degree-2 nodes
    queue<int> deg2Nodes;
    for (int i = 0; i < n; i++) {
        degree[i] = graph[i].size();
        if (degree[i] == 2) {
            deg2Nodes.push(i);
        }
    }
    
    // Process degree-2 nodes using BFS-like approach
    while (!deg2Nodes.empty()) {
        int node = deg2Nodes.front();
        deg2Nodes.pop();
        
        if (removed[node] || degree[node] != 2) continue;
        
        // Get the two neighbors
        int n1 = -1, n2 = -1, w1 = 0, w2 = 0;
        for (const auto& edge : graph[node]) {
            if (!removed[edge.to]) {
                if (n1 == -1) {
                    n1 = edge.to;
                    w1 = edge.weight;
                } else {
                    n2 = edge.to;
                    w2 = edge.weight;
                    break;
                }
            }
        }
        
        // Skip if we don't have exactly 2 valid neighbors
        if (n1 == -1 || n2 == -1 || n1 == n2) continue;
        
        int newWeight = w1 + w2;
        
        // Store the removed node's connections for path reconstruction
        replacedEdge[node] = {n1, n2};
        
        // Get existing chain between n1 and n2 (if any)
        auto minmax_pair = minmax(n1, n2);
        vector<int>& chain = edgeToChain[minmax_pair];
        chain.push_back(node);
        
        // Remove edges from neighbors to this node and update degrees
        for (auto& edge : graph[n1]) {
            if (edge.to == node) {
                edge.to = -1; // Mark for removal
                degree[n1]--;
                break;
            }
        }
        
        for (auto& edge : graph[n2]) {
            if (edge.to == node) {
                edge.to = -1; // Mark for removal
                degree[n2]--;
                break;
            }
        }
        
        // Add or update edge between n1 and n2
        bool foundEdge = false;
        for (auto& edge : graph[n1]) {
            if (edge.to == n2) {
                edge.weight = min(edge.weight, newWeight);
                foundEdge = true;
                break;
            }
        }
        if (!foundEdge) {
            graph[n1].push_back(Edge(n2, newWeight));
            degree[n1]++;
        }
        
        foundEdge = false;
        for (auto& edge : graph[n2]) {
            if (edge.to == n1) {
                edge.weight = min(edge.weight, newWeight);
                foundEdge = true;
                break;
            }
        }
        if (!foundEdge) {
            graph[n2].push_back(Edge(n1, newWeight));
            degree[n2]++;
        }
        
        // Mark node as removed
        removed[node] = true;
        degree[node] = 0;
        
        // Check if neighbors became degree-2 nodes
        if (degree[n1] == 2 && !removed[n1]) {
            deg2Nodes.push(n1);
        }
        if (degree[n2] == 2 && !removed[n2]) {
            deg2Nodes.push(n2);
        }
    }
    
    // Clean up: remove marked edges (to == -1) and clear removed nodes
    for (int i = 0; i < n; i++) {
        if (removed[i]) {
            graph[i].clear();
        } else {
            graph[i].erase(
                remove_if(graph[i].begin(), graph[i].end(),
                          [](const Edge& e) { return e.to == -1; }),
                graph[i].end()
            );
        }
    }
    
    return {graph, removed, replacedEdge, edgeToChain};
}

// Wrapper that returns just the graph (for backward compatibility)
vector<vector<Edge>> reduceGraph(const vector<vector<Edge>>& adjList) {
    return reduceGraphWithInfo(adjList).graph;
}

// Expand a path on the reduced graph back to the original graph
vector<int> expandPath(const vector<int>& reducedPath, const ReductionInfo& info) {
    if (reducedPath.size() < 2) return reducedPath;
    
    vector<int> fullPath;
    fullPath.push_back(reducedPath[0]);
    
    for (size_t i = 0; i < reducedPath.size() - 1; i++) {
        int u = reducedPath[i];
        int v = reducedPath[i + 1];
        
        // Create the edge key (always use min, max order)
        pair<int, int> edge = minmax(u, v);
        
        // Check if there were removed nodes between u and v
        auto it = info.edgeChains.find(edge);
        if (it != info.edgeChains.end()) {
            const vector<int>& chain = it->second;
            
            // Add the chain of removed nodes in correct order
            if (u < v) {
                for (int node : chain) {
                    fullPath.push_back(node);
                }
            } else {
                // Reverse order if going backwards
                for (int j = chain.size() - 1; j >= 0; j--) {
                    fullPath.push_back(chain[j]);
                }
            }
        }
        
        fullPath.push_back(v);
    }
    
    return fullPath;
}

int countActiveNodes(const vector<vector<Edge>>& graph) {
    int count = 0;
    for (const auto& adj : graph) {
        if (!adj.empty()) count++;
    }
    return count;
}

// Compact the graph by removing nodes with empty adjacency lists
// Returns pair: <compacted graph, mapping from old to new indices>
pair<vector<vector<Edge>>, vector<int>> compactGraph(const vector<vector<Edge>>& graph) {
    int n = graph.size();
    vector<int> oldToNew(n, -1);
    int newIdx = 0;
    
    // Create mapping from old to new indices
    for (int i = 0; i < n; i++) {
        if (!graph[i].empty()) {
            oldToNew[i] = newIdx++;
        }
    }
    
    // Build compacted graph
    vector<vector<Edge>> compacted(newIdx);
    for (int i = 0; i < n; i++) {
        if (!graph[i].empty()) {
            int newI = oldToNew[i];
            for (const auto& edge : graph[i]) {
                compacted[newI].push_back(Edge(oldToNew[edge.to], edge.weight));
            }
        }
    }
    
    return {compacted, oldToNew};
}

// Helper function to print the graph
void printGraph(const vector<vector<Edge>>& graph) {
    for (int i = 0; i < graph.size(); i++) {
        if (graph[i].empty()) continue;
        cout << i << ": ";
        for (const auto& edge : graph[i]) {
            cout << "(" << edge.to << ", " << edge.weight << ") ";
        }
        cout << endl;
    }
}

// BFS from a start node to find the furthest node and return the path
pair<int, vector<int>> bfsFurthest(const vector<vector<Edge>>& graph, int start) {
    int n = graph.size();
    vector<int> dist(n, -1);
    vector<int> parent(n, -1);
    queue<int> q;
    
    dist[start] = 0;
    q.push(start);
    
    int furthestNode = start;
    int maxDist = 0;
    
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        
        for (const auto& edge : graph[u]) {
            int v = edge.to;
            if (dist[v] == -1) {  // Not visited
                dist[v] = dist[u] + 1;
                parent[v] = u;
                q.push(v);
                
                if (dist[v] > maxDist) {
                    maxDist = dist[v];
                    furthestNode = v;
                }
            }
        }
    }
    
    // Reconstruct path from start to furthest node
    vector<int> path;
    int curr = furthestNode;
    while (curr != -1) {
        path.push_back(curr);
        curr = parent[curr];
    }
    reverse(path.begin(), path.end());
    
    return {furthestNode, path};
}


// Get a random active node from the graph
int getRandomActiveNode(const vector<vector<Edge>>& graph) {
    vector<int> activeNodes;
    for (int i = 0; i < graph.size(); i++) {
        if (!graph[i].empty()) {
            activeNodes.push_back(i);
        }
    }
    
    if (activeNodes.empty()) return -1;
    
    int randomIdx = rand() % activeNodes.size();
    return activeNodes[randomIdx];
}

int main() {
    srand(time(0));
    
    cout << "Loading test cases:\n";
    auto inputs = load_test_cases();
    int N, M, u, v;
    int test_case = 0;
    
    for (auto& f : inputs) {
        (*f) >> N >> M;
        
        vector<vector<Edge>> adj(N);
        vector<int> path;
        vector<int> visited(N,0);
        path.reserve(N);
        
        for (int i = 0; i < M; i++) {
            (*f) >> u >> v;
            adj[u].push_back(Edge(v, 1));
            adj[v].push_back(Edge(u, 1));
        }
        
        cout << "Processing test case " << test_case++ << ":\n";
        cout << "Graph size before reduction: " << N << " nodes, " << M << " edges\n";
        
        // Use the version with info for path reconstruction
        auto reductionInfo = reduceGraphWithInfo(adj);
        int activeNodes = countActiveNodes(reductionInfo.graph);
        int removedNodes = N - activeNodes;
        
        cout << "Graph size after reduction: " << activeNodes << " active nodes (removed " 
             << removedNodes << " nodes)\n";
        
        // BFS from random node to find furthest node
        int randomStart = getRandomActiveNode(reductionInfo.graph);
        if (randomStart != -1) {
            auto [furthestNode, reducedPath] = bfsFurthest(reductionInfo.graph, randomStart);
            
            cout << "BFS from random node " << randomStart << " to furthest node " << furthestNode << "\n";
            cout << "Path length on reduced graph: " << reducedPath.size() << " nodes\n";
            
            // Expand the path to include removed nodes
            vector<int> fullPath = expandPath(reducedPath, reductionInfo);
            cout << "Path length on original graph: " << fullPath.size() << " nodes\n";
            
         
        }
        
        // Optional: compact the graph
        auto compactResult = compactGraph(reductionInfo.graph);
        cout << "Compacted graph size: " << compactResult.first.size() << " nodes\n\n";
    }  

    return 0;
}
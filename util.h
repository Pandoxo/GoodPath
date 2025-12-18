#include <vector>
#include <queue>
#include <unordered_set>
#include <algorithm>
#include <fstream>
#include <filesystem>
#include <climits>
#include <map>
#include <random>


namespace fs=std::filesystem;
using namespace std;



struct Edge {
    int to;
    int weight;
    
    Edge(int t, int w) : to(t), weight(w) {}
};


struct ReductionInfo {
    vector<vector<Edge>> graph;           // The reduced graph
    vector<bool> removed;                 // Which nodes were removed
    vector<pair<int, int>> replacedEdge;  // For each removed node: (neighbor1, neighbor2)
    map<pair<int, int>, vector<int>> edgeChains; // Map from edge to removed nodes IN ORDER
};

// Returns the reduced graph with removal information for path reconstruction
ReductionInfo reduceGraphWithInfo(const vector<vector<Edge>>& adjList) {
    int n = adjList.size();
    vector<vector<Edge>> graph = adjList;
    vector<bool> removed(n, false);
    vector<int> degree(n);
    vector<pair<int, int>> replacedEdge(n, {-1, -1});
    
    // Map from edge (u,v) to ORDERED chain of removed nodes between them
    map<pair<int, int>, vector<int>> edgeToChain;
    
    // Initialize degrees
    for (int i = 0; i < n; i++) {
        degree[i] = graph[i].size();
    }
    
    // Process all degree-2 nodes
    bool changed = true;
    while (changed) {
        changed = false;
        
        for (int node = 0; node < n; node++) {
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
            
            // Skip if we don't have exactly 2 valid neighbors or self-loop
            if (n1 == -1 || n2 == -1 || n1 == n2) continue;
            
            int newWeight = w1 + w2;
            
            // Store the removed node's connections for path reconstruction
            replacedEdge[node] = {n1, n2};
            
            // Build the ordered chain between n1 and n2
            auto minmax_pair = minmax(n1, n2);
            vector<int> newChain;
            
            // Check if there's already a chain from n1 to node
            auto it1 = edgeToChain.find(minmax(n1, node));
            if (it1 != edgeToChain.end()) {
                // There's an existing chain between n1 and node
                if (n1 < node) {
                    newChain = it1->second;
                } else {
                    newChain = it1->second;
                    reverse(newChain.begin(), newChain.end());
                }
                edgeToChain.erase(it1);
            }
            
            // Add the current node
            newChain.push_back(node);
            
            // Check if there's already a chain from node to n2
            auto it2 = edgeToChain.find(minmax(node, n2));
            if (it2 != edgeToChain.end()) {
                // There's an existing chain between node and n2
                if (node < n2) {
                    newChain.insert(newChain.end(), it2->second.begin(), it2->second.end());
                } else {
                    auto reversed = it2->second;
                    reverse(reversed.begin(), reversed.end());
                    newChain.insert(newChain.end(), reversed.begin(), reversed.end());
                }
                edgeToChain.erase(it2);
            }
            
            // Store the new chain (always in order from min to max endpoint)
            if (n1 > n2) {
                reverse(newChain.begin(), newChain.end());
            }
            edgeToChain[minmax_pair] = newChain;
            
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
                    edge.weight = newWeight;
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
                    edge.weight = newWeight;
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
            changed = true;
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
            // Chain is stored in order from min(u,v) to max(u,v)
            if (u < v) {
                // Going forward through the chain
                for (int node : chain) {
                    fullPath.push_back(node);
                }
            } else {
                // Going backward through the chain
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




auto load_test_cases() {
    fs::path srcPath = __FILE__;
    fs::path srcDir = srcPath.parent_path();
    fs::path targetFolder = srcDir / "GP-samle-instances";

    std::vector<fs::path> filePaths;
    for (const auto& entry : fs::directory_iterator(targetFolder))
        if (entry.is_regular_file())
            filePaths.push_back(entry.path());

    std::sort(filePaths.begin(), filePaths.end(),
        [](const fs::path& a, const fs::path& b) {
            return a.filename().string() < b.filename().string();
        });

    std::vector<std::unique_ptr<std::ifstream>> inputs;

    int i = 0;
    for (const auto& p : filePaths) {
        auto f = std::make_unique<std::ifstream>(p);
        if (!f->is_open()) {
            std::cerr << "Failed to open " << p << "\n";
            continue;
        }
            inputs.push_back(std::move(f));
      
        cout<<"Loaded "<<p<<"\n";
    }
    return std::move(inputs);
}


#include <bits/stdc++.h>
using namespace std;



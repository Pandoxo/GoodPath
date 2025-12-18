#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <fstream>
#include <filesystem> 
#include <climits>
#include <queue>
#include <random>
#include <stack>
#include <unordered_set>
namespace fs = std::filesystem;

using namespace std;

struct Edge {
    int to;
    int weight;
    
    Edge(int t, int w) : to(t), weight(w) {}
};





//Greedy path . Choose next node among unmasked neighbors



vector<bool> removed(1000000,0);
vector<int> path_pool;

bool dfs_util(int node, int target,int depth,
                vector<vector<Edge>>& adj,
              vector<int>& visited,
              vector<int>& path)
{

    if (depth> 20){
        //backtrack
        path.pop_back();
        visited[node] = 0;
        for(Edge e : adj[node]){
            if(e.to != path.back())
                visited[e.to] -=1;
        }
        return false;
    }
    if (node == target) {
        return true; // found a path, stop immediately
    }

    for (Edge e : adj[node]) {
        if (visited[e.to] ==0) {
            vector<int> omit = {node,target};
            visited[e.to] =1;
            // add neighbours of our node to visited
            if (dfs_util(e.to, target, depth+1, adj, visited, path)){
                return true;  // propagate success upward
            }
        }
    }

    //backtrack
    path.pop_back();
    //unvisit neigbours - function exists on github already
    visited[node] = 0;
    for(Edge e : adj[node]){
        if(e.to != path.back())
                visited[e.to] -=1;
    }
    return false;
}


//this function assumes that neighbours of last node on a path 
bool dfs_longest_path(int node, int start,
                     vector<vector<Edge>>& adj,
                     vector<int>& visited,
                     vector<int>& path,
                     int& max_length,
                     vector<int>& best_path,int lenght = 1)
{
    path.push_back(node);
    visited[node] = 1;
    
    // Mark all neighbors of current node as visited
    for (Edge e : adj[node]) {
        visited[e.to] += 1;
    }
    
   
    
    bool found_extension = false;
    
    // Try to extend path to unvisited neighbors
    for (Edge e : adj[node]) {
        if (visited[e.to] == 1) { // Only visited once (as neighbor of current node)
            found_extension = true;
            dfs_longest_path(e.to, start, adj, visited, path, max_length, best_path,lenght + e.weight);
        }
    }

     // Update longest path if current is longer
    if (lenght > max_length) {
        max_length = lenght;
        best_path = path;
    }
    
    
    // Unmark neighbors of current node
    for (Edge e : adj[node]) {
        visited[e.to] -= 1;
    }
    
    // Unmark current node
    // Backtrack
    path.pop_back();
    visited[node]=1;
    
    return found_extension;
}

// Wrapper function to find longest path starting from a given node
int find_longest_path(int start, vector<vector<Edge>>& adj,vector<int>& visited, vector<int>& result_path) {
    int n = adj.size();
    vector<int> path;
    int max_length = 0;
    
    dfs_longest_path(start, start, adj, visited, path, max_length, result_path);
    
    return max_length;
}


vector<int> mergePathWithReroute( vector<int>& originalPath, const vector<int>& reroute) {
    if (originalPath.empty() || reroute.empty()) {
        return originalPath.empty() ? reroute : originalPath;
    }
    
    // Find where reroute starts and ends in the original path
    int rerouteStart = reroute.front();
    int rerouteEnd = reroute.back();
    
    // Find positions in original path
    int startPos = -1, endPos = -1;
    
    for (int i = 0; i < originalPath.size(); i++) {
        if (originalPath[i] == rerouteStart) startPos = i;
        if (originalPath[i] == rerouteEnd) endPos = i;
    }
    
    // Ensure startPos comes before endPos
    if (startPos > endPos) {
        swap(startPos, endPos);
    }
    
    vector<int> mergedPath;
    
    // Add path from beginning to start of reroute
    for (int i = 0; i <= startPos; i++) {
        mergedPath.push_back(originalPath[i]);
    }
    
    // Add reroute (excluding first node to avoid duplicate)
    for (int i = 1; i < reroute.size(); i++) {
        mergedPath.push_back(reroute[i]);
    }
    
    // Add remaining path from end of reroute
    for (int i = endPos + 1; i < originalPath.size(); i++) {
        mergedPath.push_back(originalPath[i]);
    }
    
    return mergedPath;
}

//find path from start to target using DFS with masking
vector<int> dfs_two_nodes(vector<vector<Edge>>& adj,int start, int target, vector<int> &visited,vector<int>& curr_path) {
    vector<int> path;
    if(curr_path.size() >= 2 ){
        //unmask neighbours of start and target
        for(Edge e : adj[start]){
            visited[e.to] -=1;
        }  
        visited[target] = 1;  
       
    }




    dfs_util(start, target,0, adj, visited, path);
    return path;  
}


// Version that returns detailed error information
struct PathValidationResult {
    bool valid;
    string error_message;
    int error_position;
};

PathValidationResult validate_path_detailed(const vector<int>& path, 
                                            const vector<vector<Edge>>& adj) {
    PathValidationResult result{true, "", -1};
    
    if (path.empty()) {
        return result;
    }
    
    if (path.size() == 1) {
        if (path[0] < 0 || path[0] >= adj.size()) {
            result.valid = false;
            result.error_message = "Node out of bounds";
            result.error_position = 0;
        }
        return result;
    }
    
    // Check for duplicate nodes
    unordered_set<int> seen;
    for (int i = 0; i < path.size(); i++) {
        if (seen.count(path[i])) {
            result.valid = false;
            result.error_message = "Duplicate node " + to_string(path[i]);
            result.error_position = i;
            return result;
        }
        seen.insert(path[i]);
    }
    
    // Check consecutive edges
    for (int i = 0; i < path.size() - 1; i++) {
        int current = path[i];
        int next = path[i + 1];
        
        if (current < 0 || current >= adj.size()) {
            result.valid = false;
            result.error_message = "Node " + to_string(current) + " out of bounds";
            result.error_position = i;
            return result;
        }
        
        if (next < 0 || next >= adj.size()) {
            result.valid = false;
            result.error_message = "Node " + to_string(next) + " out of bounds";
            result.error_position = i + 1;
            return result;
        }
        
        bool edge_exists = false;
        for (const Edge& e : adj[current]) {
            if (e.to == next) {
                edge_exists = true;
                break;
            }
        }
        
        if (!edge_exists) {
            result.valid = false;
            result.error_message = "No edge from " + to_string(current) + 
                                  " to " + to_string(next);
            result.error_position = i;
            return result;
        }
    }
    
    // Check that no non-consecutive nodes on the path are adjacent
    for (int i = 0; i < path.size(); i++) {
        int node = path[i];
        
        for (const Edge& e : adj[node]) {
            int neighbor = e.to;
            
            // Find if neighbor is on the path
            for (int j = 0; j < path.size(); j++) {
                if (path[j] == neighbor && abs(i - j) > 1) {
                    result.valid = false;
                    result.error_message = "Non-consecutive nodes " + 
                                          to_string(node) + " (position " + 
                                          to_string(i) + ") and " + 
                                          to_string(neighbor) + " (position " + 
                                          to_string(j) + ") are adjacent in graph";
                    result.error_position = i;
                    return result;
                }
            }
        }
    }
    
    return result;
}

int main(){
    ios_base::sync_with_stdio(false);
    cin.tie(0);
    srand(time(0));
    int N,M,u,v;
    cin >> N>>M;

    vector<vector<Edge>> adj(N);
    vector<int> visited(N,0);
    vector<int> path;
    path.reserve(N);
    for (int i = 0; i < M; i++) {
        cin >> u >> v;
        adj[u].push_back({v,1});
        adj[v].push_back({u,1});
    }

    vector<int> longest_path;
    srand(time(0));
    int max_lenght = find_longest_path(0, adj,visited, longest_path);

    // for(int i=0;i<10;i++){
  
    //     int start = rand()%N;
    //     fill(visited.begin(),visited.end(),0);
    //     int length = find_longest_path(start, adj,visited, longest_path);
    //     if(length > max_lenght){
    //         max_lenght = length;
    //     }
    // }
    cout << longest_path.size() << "\n";
    for(int i = 0; i< (int)longest_path.size();i++){
        if(i) cout << ' ';
        cout << longest_path[i];
    }
          
    return 0;
}
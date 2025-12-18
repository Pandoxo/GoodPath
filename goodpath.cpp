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
#include "util.h"
namespace fs = std::filesystem;

using namespace std;





void sort_by_bfs(vector<vector<Edge>>& adj, int N) {
    queue<int> q;
    vector<bool> visited(N, false);
    vector<int> distances(N,0);
    vector<int> parents(N,-1);
    q.push(0);
    visited[0] = true;

    while(!q.empty()) {
        int node = q.front();
        q.pop();

        for(Edge e : adj[node]) {
            if(!visited[e.to]) {
                visited[e.to] = true;
                //parents[neighbour] = node;
                distances[e.to] +=1;
                q.push(e.to);
            }
        }
    }
    //reconstruct path
    // vector<int> path;
    // int current = N-1;
    // while(current != -1) {
    //     path.push_back(current);
    //     current = parents[current];
    // }
    // reverse(path.begin(), path.end());
    vector<vector<Edge>> temp = adj;
    for (vector<Edge> &neighbours : adj)
    {
        sort(neighbours.begin(), neighbours.end(), [&](const Edge &a, const Edge &b)
             { return temp[a.to].size() < temp[b.to].size(); });
    }
  
    
}




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

int i=0;
//this function assumes that neighbours of last node on a path 
void dfs_longest_path(int& node,
                     vector<vector<Edge>>& adj,
                     vector<int>& visited,
                     vector<int>& path,
                     int& max_length,
                     vector<int>& best_path,
                     chrono::steady_clock::time_point& start_time,
                     int search_duration,
                     int lenght = 1)
{
    // Check if time limit exceeded
    auto now = chrono::steady_clock::now();
    if (chrono::duration_cast<chrono::seconds>(now - start_time).count() > search_duration) {
        return;
    }

    path.push_back(node);
    visited[node] = 1;

    
    // Mark all neighbors of current node as visited
    for (Edge e : adj[node]) {
        visited[e.to] += 1;
    }
    
    
    // Try to extend path to unvisited neighbors
    for (Edge e : adj[node]) {
        if (visited[e.to] == 1) { // Only visited once (as neighbor of current node)
            dfs_longest_path(e.to, adj, visited, path, max_length, best_path, start_time, search_duration, lenght + e.weight);
        }
    }

     // Update longest path if current is longer
    if (lenght > max_length) {
        max_length = lenght;
        best_path = path;
    }
   
    
    // Backtrack
    // Unmark neighbors of current node
    for (Edge e : adj[node]) {
        visited[e.to] -= 1;
    }
    
    // Unmark current node
    path.pop_back();
    
    visited[node] = 1;
    
}

// Wrapper function to find longest path starting from a given node
vector<int> find_longest_path(int start, vector<vector<Edge>>& adj, vector<int>& visited ) {
    vector<int> path;
    int max_length = 0;
    vector<int> result_path;
    
    auto start_time = chrono::steady_clock::now();
    int search_duration = 20; 
    
    dfs_longest_path(start, adj, visited, path, max_length, result_path, start_time, search_duration);
    
   return result_path;
}
    
void update_visited(vector<vector<Edge>>& adj,vector<int>& path,vector<int>& visited){
    fill(visited.begin(), visited.end(), 0);
    
    for(int i=1; i< (int)path.size();i++){
        int node = path[i];
        for(Edge e: adj[node]){
            visited[e.to] =+1;
        }
    }
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

    srand(time(0));

    cout <<"Loading test cases:\n";
    auto inputs = load_test_cases();
    int N,M,u,v;
    int test_case = -2;
    
    for (auto& f : inputs) {
        (*f) >> N >> M;

        vector<vector<Edge>> adj(N);
        vector<int> visited(N,0);
        vector<int> path;
        path.reserve(N);
        
        for (int i = 0; i < M; i++) {
            (*f) >> u >> v;
            adj[u].push_back({v,1});
            adj[v].push_back({u,1});
        }

        vector<int> longest_path;
        vector<int> from_zero_path;

        srand(time(0));
        //auto reduced = reduceGraphWithInfo(adj);
        int start = rand()%N;
        cout << "Test case: " << test_case++ <<  " size: "<< adj.size()<< "\n";
        sort_by_bfs(adj,N);
    
        longest_path = find_longest_path(start, adj,visited);
        update_visited(adj,longest_path,visited);
        from_zero_path = find_longest_path(start, adj,visited);   

        cout << "Longest path length from rand: " << longest_path.size() << endl;
        cout << "extend tail: " << from_zero_path.size() << endl;


        // for(int i=0;i<10;i++){
        //     //int start = getRandomActiveNode(adj);
        //     int start = rand()%N;
        //     fill(visited.begin(),visited.end(),0);
        //     int length = find_longest_path(start, adj,visited, longest_path);
        //     if(length > max_lenght){
        //         max_lenght = length;
        //     }
        // }
        //path after expansion
       // path = expandPath(longest_path,reduced);
        //cout << "Longest path length from rand: " << longest_path.size() << endl;
        
        PathValidationResult result = validate_path_detailed(longest_path,adj);
        cout << result.error_message << "\n";
        
    }  
    return 0;
}
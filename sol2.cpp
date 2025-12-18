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
                     chrono::steady_clock::time_point& last_print_time,
                     int lenght = 1)
{
    // Check if time limit exceeded
    auto now = chrono::steady_clock::now();
    float elapsed_seconds = chrono::duration<float>(now - start_time).count();

    if (elapsed_seconds > search_duration) {
        return;
    }

    // Print max_length every 0.5 seconds
    // if (chrono::duration_cast<chrono::milliseconds>(now - last_print_time).count() >= 500) {
    //     cout << "Time: " << elapsed_seconds << " max_length: " << max_length << endl;
    //     last_print_time = now;
    // }

    path.push_back(node);
    visited[node] = 1;

    
    // Mark all neighbors of current node as visited
    for (Edge e : adj[node]) {
        visited[e.to] += 1;
    }
    
    
    // Try to extend path to unvisited neighbors
    for (Edge e : adj[node]) {
        if (visited[e.to] == 1) { // Only visited once (as neighbor of current node)
            dfs_longest_path(e.to, adj, visited, path, max_length, best_path, start_time, search_duration, last_print_time, lenght + e.weight);
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
vector<int> find_longest_path(int start, vector<vector<Edge>>& adj, vector<int>& visited,int search_duration  ) {
    vector<int> path;
    int max_length = 0;
    vector<int> result_path;
    
    auto start_time = chrono::steady_clock::now();
    auto last_print = chrono::steady_clock::now();
    dfs_longest_path(start, adj, visited, path, max_length, result_path, start_time, search_duration,last_print);
    
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

vector<int> merge_paths_from_same_start(const vector<int>& path1, const vector<int>& path2) {
    if (path1.empty()) return path2;
    if (path2.empty()) return path1;
    
    vector<int> result;
    
    // Determine which path is longer
    const vector<int>& longer = path1.size() >= path2.size() ? path1 : path2;
    const vector<int>& shorter = path1.size() >= path2.size() ? path2 : path1;
    
    // Reverse the shorter path (excluding the start node)
    for (int i = shorter.size() - 1; i >= 1; i--) {
        result.push_back(shorter[i]);
    }
    
    // Add the longer path (including the start node)
    for (int i = 0; i < longer.size(); i++) {
        result.push_back(longer[i]);
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
    for (int i = 0; i < M; i++) {
        cin >> u >> v;
        adj[u].push_back({v,1});
        adj[v].push_back({u,1});
    }

    vector<int> path;
    vector<int> longest_path;
    
    srand(time(0));
    //auto reduced = reduceGraphWithInfo(adj);
    int start = rand()%N;
    int longest_start;
    sort_by_bfs(adj,N);
    
    for(int i=0;i<10;i++){
        //int start = getRandomActiveNode(adj);
        start = rand()%N;
        fill(visited.begin(),visited.end(),0);
        path = find_longest_path(start, adj,visited,1);
        if(path.size() > longest_path.size()){
            longest_path = path;
            longest_start = start;
        }
    }
    vector<int> from_zero_path;
    update_visited(adj,longest_path,visited);
    from_zero_path = find_longest_path(longest_start, adj,visited,1);   
    vector<int> result = merge_paths_from_same_start(longest_path,from_zero_path);

    cout << result.size() << "\n";
    for(int i = 0; i< (int)result.size();i++){
        if(i) cout << ' ';
        cout << result[i];
    }
          
    return 0;
}
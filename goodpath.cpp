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


struct Edge{
    int to;
    int len;
    int path_start;
    int path_len;
};



vector<int> bfs(vector<vector<int>>& adj, int N) {
    queue<int> q;
    vector<bool> visited(N, false);
    //vector<int> distances(N,0);
    vector<int> parents(N,-1);
    q.push(0);
    visited[0] = true;
    int count = 1;  

    while(!q.empty()) {
        int node = q.front();
        q.pop();

        for(int neighbour : adj[node]) {
            if(!visited[neighbour]) {
                visited[neighbour] = true;
                parents[neighbour] = node;
                q.push(neighbour);
                count++;
            }
        }
    }
    //reconstruct path
    vector<int> path;
    int current = N-1;
    while(current != -1) {
        path.push_back(current);
        current = parents[current];
    }
    reverse(path.begin(), path.end());
    return path;
    
}


//Greedy path . Choose next node among unmasked neighbors
int extend_path(
    vector<vector<int>> &adj,
    vector<int>& path,
    unordered_set<int>& mask) {

    int current_node = path.back();

    int next_node = -1;
    int min_degree = INT_MAX;

    for(int neighbour : adj[current_node]) {
        if(mask.find(neighbour) == mask.end()) {
            next_node = neighbour;
            break;
        }
    }

    if(next_node != -1) {
        for(int neighbour : adj[current_node]) {
            mask.insert(neighbour);
        }
        path.push_back(next_node);
    }else{
        mask.insert(current_node);
        path.pop_back();
    }
    return next_node;
}


vector<bool> removed(1000000,0);
vector<int> path_pool;


void reduce_degree_two(vector<vector<Edge>> &adj,int n ) {
    queue<int> q;

    for (int v = 0; v < n; v++) {
        if (!removed[v] && adj[v].size() == 2)
            q.push(v);
    }

    while (!q.empty()) {
        int v = q.front(); q.pop();
        if (removed[v] || adj[v].size() != 2) continue;

        Edge e1 = adj[v][0];
        Edge e2 = adj[v][1];
        int u = e1.to;
        int w = e2.to;

        removed[v] = true;

        // Remove v from neighbors
        auto remove_v = [&](int x){
            adj[x].erase(remove_if(adj[x].begin(), adj[x].end(),
                    [&](Edge &e){ return e.to == v; }), adj[x].end());
        };
        remove_v(u);
        remove_v(w);

        if (u != w) { // avoid self-loop
            // Store new edge in global path pool
            int start_idx = path_pool.size();

            // Append paths from e1
            if (e1.path_len > 0)
                path_pool.insert(path_pool.end(), path_pool.begin() + e1.path_start,
                                 path_pool.begin() + e1.path_start + e1.path_len);

            // Add the removed node itself
            path_pool.push_back(v);

            // Append paths from e2
            if (e2.path_len > 0)
                path_pool.insert(path_pool.end(), path_pool.begin() + e2.path_start,
                                 path_pool.begin() + e2.path_start + e2.path_len);

            int total_len = e1.len + e2.len;
            int total_path_len = path_pool.size() - start_idx;

            // Add edge u -> w
            adj[u].push_back({w, total_len, start_idx, total_path_len});
            adj[w].push_back({u, total_len, start_idx, total_path_len});
        }

        // Checout << "Test case "<<test_case++<<": N = "<<N<<", M = "<<M<<", deg2 = "<<count<<"\n";ck neighbors
        if (!removed[u] && adj[u].size() == 2) q.push(u);
        if (!removed[w] && adj[w].size() == 2) q.push(w);

        adj[v].clear();
    }
}


// vector<int> find_cycle(
//     vector<vector<int>> &adj,
//     int start_node) {

//     stack<pair<int,vector<int>>> s;
//     s.push({start_node,{}});
//     vector<int> path;

//     while (!s.empty()) {
//         auto [node,path] = s.top();
//         s.pop();
//         for(int neighbour : adj[node]){
//             if(!check_if_neighbor_in_path(adj,path,neighbour)){
//                 path.push_back(node);
//                 path.push_back(neighbour);
//                 return path;
//         }
//     }
//     return {};
// }
// // Reconstruct original path from reduced path of nodes
// vector<int> reconstruct_path( vector<vector<Edge>> &adj,const vector<int> &reduced_path) {
//     vector<int> full_path;
//     for (size_t i = 0; i + 1 < reduced_path.size(); i++) {
//         int u = reduced_path[i];
//         int v = reduced_path[i+1];

//         auto it = find_if(adj[u].begin(), adj[u].end(), [&]( Edge &e){ return e.to == v; });
//         if (it != adj[u].end()) {
//             full_path.push_back(u);
//             if (it->path_len > 0)
//                 full_path.insert(full_path.end(),
//                     path_pool.begin() + it->path_start,
//                     path_pool.begin() + it->path_start + it->path_len);
//         } else {
//             full_path.push_back(u);
//         }
//     }
//     full_path.push_back(reduced_path.back());
//     return full_path;
// }

bool dfs_util(int node, int target,int depth,
                vector<vector<int>>& adj,
              unordered_set<int> mask,
              vector<int>& path)
{
    mask.insert(node);
    path.push_back(node);

    if (depth> 20){
        //backtrack
        path.pop_back();
        mask.erase(node);
        return false;
    }
    if (node == target) {
        return true; // found a path, stop immediately
    }

    for (int nei : adj[node]) {
        if (mask.find(nei) == mask.end()) {
            vector<int> omit = {node,target};
            mask_neighbours(adj,node,mask,omit);
            if (dfs_util(nei, target, depth+1, adj, mask, path)){
                return true;  // propagate success upward
            }
        }
    }

    //backtrack
    path.pop_back();
    mask.erase(node); 
    return false;
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
vector<int> dfs_two_nodes(vector<vector<int>>& adj,int start, int target, unordered_set<int>& mask,vector<int>& curr_path) {
    vector<int> path;
    int omit_index = -1;
    if(curr_path.size() >= 2 ){
        int start_i_on_path = find(curr_path.begin(), curr_path.end(), start) - curr_path.begin(); 
        omit_index= curr_path[start_i_on_path - 1];
        unmask_neighbours(adj,mask,start,omit_index);
        
        int target_i_on_path = find(curr_path.begin(), curr_path.end(), target) - curr_path.begin();
        omit_index = curr_path[max(target_i_on_path + 1,(int)path.size()-1)];
        unmask_neighbours(adj,mask,target,omit_index);
        mask.erase(start);
        mask.erase(target);
    }

      cout << "\nMask after:\n";
        for(int m : mask) {
            cout << m << " ";
        }



    dfs_util(start, target,0, adj, mask, path);
    return path;  
}

int main(){

    srand(time(0));

    cout <<"Loading test cases:\n";
    auto inputs = load_test_cases();
    int N,M,u,v;
    int test_case = 0;
    
    for (auto& f : inputs) {
        (*f) >> N >> M;

        vector<vector<int>> adj(N);
        unordered_set<int> mask;
        vector<int> path;
        path.reserve(N);
        
        for (int i = 0; i < M; i++) {
            (*f) >> u >> v;
            adj[u].push_back(v);
            adj[v].push_back(u);
        }
        
        path = bfs(adj,N);
        mask_neighbours_on_path(adj,path,mask,0,path.size());
        cout<<"BFS Path: ";
        for(int node : path) {
            cout << node << " ";
        }
         cout << "\nMask before:\n";
        for(int m : mask) {
            cout << m << " ";
        }
       
        unmask_neighbours_on_path(adj,path,mask,2,6);
        for(int i = 2; i<6;i++){
            mask.erase(path[i]);
        }
         cout << "\nMask after:\n";
        for(int m : mask) {
            cout << m << " ";
        }

        cout << "\nDFS Path: ";

        vector<int> curr_path = dfs_two_nodes(adj,2,6,mask,path);
        for(int node : curr_path) {
            cout << node << " ";
        }  
        cout << "\n";
        break;
        //cout<<"Test case "<<test_case++<<": Path length = "<<path.size()<<"\n";
    }  
    return 0;
}
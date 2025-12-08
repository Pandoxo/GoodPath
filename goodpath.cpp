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


struct Edge{
    int to;
    int len;
    int path_start;
    int path_len;
};

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
    return path;
    
}

// mask neighbours of nodes in path from start_index to end_index
void mask_neighbours(
    vector<vector<int>> &adj,
    vector<int>& path,
    unordered_set<int>& mask,
    int start_index = 0,
    int end_index = -1 ) {

    for(int i=start_index; i < end_index;i++){
        int node = path[i];
        for(int neighbour : adj[node]) {
            if(mask.find(neighbour) == mask.end()) {
                mask.insert(neighbour);
            }
        }
    }
}

//unamsk neighbors of a node, omitting certain index
void unmask_neighbours(
    vector<vector<int>> &adj,
    vector<int>& path,
    unordered_set<int>& mask,
    int node,
    int omit_index = -1) {

    for(int neighbour : adj[node]) {
        if(mask.find(neighbour) != mask.end() &&
            neighbour != omit_index) {
            mask.erase(neighbour);
        }
    }
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

// Reconstruct original path from reduced path of nodes
vector<int> reconstruct_path( vector<vector<Edge>> &adj,const vector<int> &reduced_path) {
    vector<int> full_path;
    for (size_t i = 0; i + 1 < reduced_path.size(); i++) {
        int u = reduced_path[i];
        int v = reduced_path[i+1];

        auto it = find_if(adj[u].begin(), adj[u].end(), [&]( Edge &e){ return e.to == v; });
        if (it != adj[u].end()) {
            full_path.push_back(u);
            if (it->path_len > 0)
                full_path.insert(full_path.end(),
                    path_pool.begin() + it->path_start,
                    path_pool.begin() + it->path_start + it->path_len);
        } else {
            full_path.push_back(u);
        }
    }
    full_path.push_back(reduced_path.back());
    return full_path;
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
    
        cout << "Test case " << test_case ++<< " nodes: " << N 
            << ", edges: " << M << "\n";


      
        
        //cout<<"Test case "<<test_case++<<": Path length = "<<path.size()<<"\n";
    
    }
    
    return 0;
}
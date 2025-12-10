#include <vector>
#include <queue>
#include <unordered_set>
#include <algorithm>
#include <fstream>
#include <filesystem>
#include <climits>
#include <random>


namespace fs=std::filesystem;
using namespace std;



//check if any neighbor of node is in path
bool check_if_neighbor_in_path(
    vector<vector<int>> &adj,
    vector<int>& path,
    int node) {

    for(int neighbour : adj[node]) {
        if(find(path.begin(), path.end(), neighbour) != path.end()) {
            return true;
        }
    }
    return false;
}

//check if any neighbor of node is in mask
bool check_if_neighbor_in_mask(
    vector<vector<int>> &adj,
    unordered_set<int>& mask,
    int node) {

    for(int neighbour : adj[node]) {
        if(mask.find(neighbour) != mask.end()) {
            return true;
        }
    }
    return false;
}


//unamsk neighbors of a node, omitting certain index
void unmask_neighbours(
    vector<vector<int>> &adj,
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
// unmask neighbours of nodes in path from start_index to end_index (exclusive)
void unmask_neighbours_on_path(
    vector<vector<int>> &adj,
    vector<int>& path,
    unordered_set<int>& mask,
    int start_index = 0,
    int end_index = -1 ) {

    for(int i=start_index; i < end_index;i++){
        int node = path[i];
        for(int neighbour : adj[node]) {
            if(mask.find(neighbour) == mask.end() &&
             find(path.begin(), path.end(), neighbour) == path.end()) {
                mask.erase(neighbour);
            }
        }
    }
}

// mask neighbours of nodes in path from start_index to end_index (exclusive)
void mask_neighbours_on_path(
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
// mask neighbours of nodes in path from start_index to end_index (exclusive)
void mask_neighbours(
    vector<vector<int>> &adj,
    int node,
    unordered_set<int>& mask,
    vector<int>& omit) {

    for(int neighbour : adj[node]) {
        if(mask.find(neighbour) == mask.end() && find(omit.begin(), omit.end(), neighbour) == omit.end()) {
            mask.insert(neighbour);
        }
    }
    
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



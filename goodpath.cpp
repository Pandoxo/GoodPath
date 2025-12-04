#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <fstream>
#include <filesystem> 
namespace fs = std::filesystem;

using namespace std;


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

int main(){

    cout <<"Loading test cases:\n";
    auto inputs = load_test_cases();
    int N,M,u,v;
    *inputs[0] >> N>>M;
    int test_case = 0;
    for (auto& f : inputs) {
        vector<vector<int>> adj(N);
        cout << "Test case " << test_case++ << ":\n";
        (*f) >> N >> M;
        for (int i = 0; i < M; i++) {
            (*f) >> u >> v;
            adj[u].push_back(v);
            adj[v].push_back(u);

        }
    }
    
    return 0;
}
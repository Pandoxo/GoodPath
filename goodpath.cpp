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

void sort_by_bfs(vector<vector<Edge>> &adj, int start, int N)
{
    queue<int> q;
    vector<bool> visited(N, false);
    vector<int> distances(N, 0);
    vector<int> parents(N, -1);
    q.push(start);
    visited[start] = true;

    while (!q.empty())
    {
        int node = q.front();
        q.pop();

        for (Edge e : adj[node])
        {
            if (!visited[e.to])
            {
                visited[e.to] = true;
                // parents[neighbour] = node;
                distances[e.to] += 1;
                q.push(e.to);
            }
        }
    }

    vector<vector<Edge>> temp = adj;

    for (vector<Edge> &neighbours : adj)
    {
        sort(neighbours.begin(), neighbours.end(),
             [&](const Edge &a, const Edge &b)
             {
                 if (temp[a.to].size() == temp[b.to].size())
                 {
                     return distances[a.to] < distances[b.to];
                 }
                 return temp[a.to].size() < temp[b.to].size();
             });
    }
}

// this function assumes that neighbours of last node on a path
void dfs_longest_path(int &node,
                      vector<vector<Edge>> &adj,
                      vector<int> &visited,
                      vector<int> &path,
                      int &max_length,
                      vector<int> &best_path,
                      chrono::steady_clock::time_point &start_time,
                      int search_duration,
                      chrono::steady_clock::time_point &last_print_time,
                      int lenght = 1)
{
    // Check if time limit exceeded
    auto now = chrono::steady_clock::now();
    float elapsed_seconds = chrono::duration<float>(now - start_time).count();

    if (elapsed_seconds > search_duration)
    {
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
    for (Edge e : adj[node])
    {
        visited[e.to] += 1;
    }

    // Try to extend path to unvisited neighbors
    for (Edge e : adj[node])
    {
        if (visited[e.to] == 1)
        { // Only visited once (as neighbor of current node)
            dfs_longest_path(e.to, adj, visited, path, max_length, best_path, start_time, search_duration, last_print_time, lenght + e.weight);
        }
    }

    // Update longest path if current is longer
    if (lenght > max_length)
    {
        max_length = lenght;
        best_path = path;
    }

    // Backtrack
    // Unmark neighbors of current node
    for (Edge e : adj[node])
    {
        visited[e.to] -= 1;
    }

    // Unmark current node
    path.pop_back();

    visited[node] = 1;
}

// Wrapper function to find longest path starting from a given node
vector<int> find_longest_path(int start, vector<vector<Edge>> &adj, vector<int> &visited)
{
    vector<int> path;
    int max_length = 0;
    vector<int> result_path;

    auto start_time = chrono::steady_clock::now();
    int search_duration = 1;
    auto last_print = chrono::steady_clock::now();
    dfs_longest_path(start, adj, visited, path, max_length, result_path, start_time, search_duration, last_print);

    return result_path;
}

void dfs_blind_util(int &node,
                    vector<vector<Edge>> &adj,
                    vector<int> &visited,
                    vector<bool> &was_in_path,
                    vector<int> &path,
                    int &max_length,
                    vector<int> &best_path,
                    int lenght = 1)
{
    // Check if time limit exceeded

    path.push_back(node);
    visited[node] = 1;
    was_in_path[node] = true;

    // Mark all neighbors of current node as visited
    for (Edge e : adj[node])
    {
        visited[e.to] += 1;
    }

    // Try to extend path to unvisited neighbors
    for (Edge e : adj[node])
    {
        if (visited[e.to] == 1 && !was_in_path[e.to])
        { // Only visited once (as neighbor of current node)
            dfs_blind_util(e.to, adj, visited, was_in_path, path, max_length, best_path, lenght + e.weight);
        }
    }

    // Update longest path if current is longer
    if (lenght > max_length)
    {
        max_length = lenght;
        best_path = path;
    }

    // Backtrack
    // Unmark neighbors of current node
    for (Edge e : adj[node])
    {
        visited[e.to] -= 1;
    }

    // Unmark current node
    path.pop_back();

    visited[node] = 1;
}

// Wrapper function to find longest path starting from a given node
vector<int> dfs_blind(int start, vector<vector<Edge>> &adj, vector<int> &visited)
{
    int max_length = 0;
    vector<int> result_path;
    vector<int> path;
    vector<bool> was_in_path(adj.size(), false);
    dfs_blind_util(start, adj, visited, was_in_path, path, max_length, result_path);

    return result_path;
}

vector<int> dfs_blind_extend(int start, vector<vector<Edge>> &adj, vector<int> &path, vector<int> &visited)
{
    int max_length = 0;
    vector<int> result_path;
    vector<int> new_path;
    vector<bool> was_in_path(adj.size(), false);
    for (int node : path)
        was_in_path[node] = true;

    dfs_blind_util(start, adj, visited, was_in_path, new_path, max_length, result_path);

    return result_path;
}

std::vector<int> trimTail(vector<vector<Edge>> adj, const std::vector<int> &vec, vector<int> &visited, double percent)
{
    // Validate input
    if (percent < 0.0 || percent > 100.0)
    {
        throw std::invalid_argument("Percent must be between 0 and 100");
    }

    if (vec.empty())
    {
        return vec;
    }

    // Calculate how many elements to keep
    size_t numToKeep = static_cast<size_t>(
        std::ceil(vec.size() * (100.0 - percent) / 100.0));
    // unmark nodes in the trimmed tail
    for (int i = numToKeep; i < vec.size(); i++)
    {
        int node = vec[i];
        visited[node] -= 1;
        for (Edge e : adj[node])
        {
            visited[e.to] -= 1;
        }
    }
    // Ensure we keep at least 0 elements and at most all elements
    numToKeep = std::min(numToKeep, vec.size());

    // Create new vector with only the first numToKeep elements
    return std::vector<int>(vec.begin(), vec.begin() + numToKeep);
}

void dfs_reverse_tail(int &node,
                      vector<vector<Edge>> &adj,
                      vector<int> &visited,
                      vector<int> &path,
                      vector<int> &result,
                      int lenght, bool &found)
{

    // Print max_length every 0.5 seconds
    // if (chrono::duration_cast<chrono::milliseconds>(now - last_print_time).count() >= 500) {
    //     cout << "Time: " << elapsed_seconds << " max_length: " << max_length << endl;
    //     last_print_time = now;
    // }

    if (found)
    {
        return;
    }

    path.push_back(node);
    visited[node] = 1;

    // Mark all neighbors of current node as visited
    for (Edge e : adj[node])
    {
        visited[e.to] += 1;
    }

    // Try to extend path to unvisited neighbors
    for (Edge e : adj[node])
    {
        if (visited[e.to] > 1 && find(path.begin(), path.end(), e.to) == path.end() && lenght > 1000)
        {
            path.push_back(e.to);
            result = path;
            found = true;
            return;
        }
        else if (visited[e.to] == 1)
        { // Only visited once (as neighbor of current node)
            dfs_reverse_tail(e.to, adj, visited, path, result, lenght + e.weight, found);
        }
    }

    // Update longest path if current is longer

    // Backtrack
    // Unmark neighbors of current node
    for (Edge e : adj[node])
    {
        visited[e.to] -= 1;
    }

    // Unmark current node
    path.pop_back();

    visited[node] = 1;
}

void update_visited(vector<vector<Edge>> &adj, vector<int> &path, vector<int> &visited)
{

    fill(visited.begin(), visited.end(), 0);

    for (int i = 1; i < (int)path.size(); i++)
    {
        int node = path[i];
        visited[node] = 1;
        for (Edge e : adj[node])
        {
            visited[e.to] += 1;
        }
    }
}

vector<int> reverseTail(vector<vector<Edge>> &adj, vector<int> &path, vector<int> &visited)
{
    vector<int> trimmed = trimTail(adj, path, visited, 10.0);
    cout << "path_size after trim " << trimmed.size() << "\n";
    if (trimmed.empty())
    {
        return path;
    }
    fill(visited.begin(), visited.end(), 0);
    update_visited(adj, trimmed, visited);

    int start = trimmed.back();
    for (Edge e : adj[start])
    {
        visited[e.to] -= 1;
    }
    trimmed.pop_back();

    vector<int> result;
    bool found = false;

    // Find tail extension starting from the last node of trimmed path
    dfs_reverse_tail(start, adj, visited, trimmed, result, 0, found);
    if (result.size() == 0)
        return {};
    int node_before_turning_point = result.back();
    trimmed = result;
    int turning_point = -1;
    vector<int>::iterator pos;
    for (Edge e : adj[node_before_turning_point])
    {
        pos = find(trimmed.begin(), trimmed.end(), e.to);
        if (pos != trimmed.end() && pos != trimmed.end() - 2)
        {
            turning_point = *pos;
            break;
        }
    }

    if (turning_point == -1)
        return {};

    reverse(trimmed.begin() + turning_point + 1, trimmed.end());
    trimmed.pop_back();

    fill(visited.begin(), visited.end(), 0);
    update_visited(adj, trimmed, visited);
    return trimmed;
}

vector<int> merge_paths_from_same_start(const vector<int> &path1, const vector<int> &path2)
{
    if (path1.empty())
        return path2;
    if (path2.empty())
        return path1;

    vector<int> result;

    // Determine which path is longer
    const vector<int> &longer = path1.size() >= path2.size() ? path1 : path2;
    const vector<int> &shorter = path1.size() >= path2.size() ? path2 : path1;

    // Reverse the shorter path (excluding the start node)
    for (int i = shorter.size() - 1; i >= 1; i--)
    {
        result.push_back(shorter[i]);
    }

    // Add the longer path (including the start node)
    for (int i = 0; i < longer.size(); i++)
    {
        result.push_back(longer[i]);
    }

    return result;
}

vector<int> mergePathWithReroute(vector<int> &originalPath, const vector<int> &reroute)
{
    if (originalPath.empty() || reroute.empty())
    {
        return originalPath.empty() ? reroute : originalPath;
    }

    // Find where reroute starts and ends in the original path
    int rerouteStart = reroute.front();
    int rerouteEnd = reroute.back();

    // Find positions in original path
    int startPos = -1, endPos = -1;

    for (int i = 0; i < originalPath.size(); i++)
    {
        if (originalPath[i] == rerouteStart)
            startPos = i;
        if (originalPath[i] == rerouteEnd)
            endPos = i;
    }

    // Ensure startPos comes before endPos
    if (startPos > endPos)
    {
        swap(startPos, endPos);
    }

    vector<int> mergedPath;

    // Add path from beginning to start of reroute
    for (int i = 0; i <= startPos; i++)
    {
        mergedPath.push_back(originalPath[i]);
    }

    // Add reroute (excluding first node to avoid duplicate)
    for (int i = 1; i < reroute.size(); i++)
    {
        mergedPath.push_back(reroute[i]);
    }

    // Add remaining path from end of reroute
    for (int i = endPos + 1; i < originalPath.size(); i++)
    {
        mergedPath.push_back(originalPath[i]);
    }

    return mergedPath;
}

// Version that returns detailed error information
struct PathValidationResult
{
    bool valid;
    string error_message;
    int error_position;
};

PathValidationResult validate_path_detailed(const vector<int> &path,
                                            const vector<vector<Edge>> &adj)
{
    PathValidationResult result{true, "", -1};

    if (path.empty())
    {
        return result;
    }

    if (path.size() == 1)
    {
        if (path[0] < 0 || path[0] >= adj.size())
        {
            result.valid = false;
            result.error_message = "Node out of bounds";
            result.error_position = 0;
        }
        return result;
    }

    // Check for duplicate nodes
    unordered_set<int> seen;
    for (int i = 0; i < path.size(); i++)
    {
        if (seen.count(path[i]))
        {
            result.valid = false;
            result.error_message = "Duplicate node " + to_string(path[i]);
            result.error_position = i;
            return result;
        }
        seen.insert(path[i]);
    }

    // Check consecutive edges
    for (int i = 0; i < path.size() - 1; i++)
    {
        int current = path[i];
        int next = path[i + 1];

        if (current < 0 || current >= adj.size())
        {
            result.valid = false;
            result.error_message = "Node " + to_string(current) + " out of bounds";
            result.error_position = i;
            return result;
        }

        if (next < 0 || next >= adj.size())
        {
            result.valid = false;
            result.error_message = "Node " + to_string(next) + " out of bounds";
            result.error_position = i + 1;
            return result;
        }

        bool edge_exists = false;
        for (const Edge &e : adj[current])
        {
            if (e.to == next)
            {
                edge_exists = true;
                break;
            }
        }

        if (!edge_exists)
        {
            result.valid = false;
            result.error_message = "No edge from " + to_string(current) +
                                   " to " + to_string(next);
            result.error_position = i;
            return result;
        }
    }

    // Check that no non-consecutive nodes on the path are adjacent
    for (int i = 0; i < path.size(); i++)
    {
        int node = path[i];

        for (const Edge &e : adj[node])
        {
            int neighbor = e.to;

            // Find if neighbor is on the path
            for (int j = 0; j < path.size(); j++)
            {
                if (path[j] == neighbor && abs(i - j) > 1)
                {
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

// dfs from the end of the path

int countDuplicateOccurrences(const std::vector<int> &vec)
{
    std::unordered_map<int, int> freq;

    for (int num : vec)
    {
        freq[num]++;
    }

    int duplicateCount = 0;
    for (const auto &pair : freq)
    {
        if (pair.second > 1)
        {
            duplicateCount += pair.second - 1; // Count extras beyond first
        }
    }

    return duplicateCount;
}

void shuffle_adj(vector<vector<Edge>> &adj)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    for (vector<Edge> &neighbours : adj)
    {
        shuffle(neighbours.begin(), neighbours.end(), gen);
    }
}

int main()
{

    srand(164205);

    cout << "Loading test cases:\n";
    auto inputs = load_test_cases();
    int N, M, u, v;
    int test_case = 1;

    for (auto &f : inputs)
    {
        (*f) >> N >> M;

        vector<vector<Edge>> adj(N);
        vector<int> visited(N, 0);

        for (int i = 0; i < M; i++)
        {
            (*f) >> u >> v;
            adj[u].push_back({v, 1});
            adj[v].push_back({u, 1});
        }

        vector<int> longest_path;
        vector<int> longest_path_blind;

        vector<int> path;
        vector<int> from_zero_path;
        vector<int> result;
        vector<vector<int>> solutions;

        int start;
        srand(time(0));
        // auto reduced = reduceGraphWithInfo(adj);
        // int start = getRandomActiveNode(reduced.graph);
        // sort_by_bfs(reduced.graph,start,N);
        cout << "Test case: " << test_case++ << '\n';
        auto compare_by_size =
            ([&](const vector<int> a, const vector<int> &b)

             { return a.size() > b.size(); });

        path = find_longest_path(0, adj, visited);
        cout << path.size() << "\n";
        path = reverseTail(adj, path, visited);
        fill(visited.begin(), visited.end(), 0);
        update_visited(adj, path, visited);
        from_zero_path = dfs_blind_extend(path.back(), adj, path, visited);
        path = merge_paths_from_same_start(path, from_zero_path);
        cout << path.size() << "\n";

        // for (int i = 0; i < 10; i++)
        // {

        //     shuffle_adj(adj);
        //     path = dfs_blind(start, adj, visited);
        //     fill(visited.begin(), visited.end(), 0);
        //     update_visited(adj, path, visited);

        //     from_zero_path = dfs_blind_extend(start, adj, path, visited);
        //     path = merge_paths_from_same_start(path, from_zero_path);
        //     solutions.push_back(path);
        //     sort(solutions.begin(), solutions.end(), compare_by_size);
        //     if (solutions.size() > 5)
        //     {
        //         solutions.pop_back();
        //     }
        // }
        // fill(visited.begin(), visited.end(), 0);
        // for (int i = 0; i < 8; i++)
        // {
        //     // int start = getRandomActiveNode(adj);
        //     start = rand() % N;
        //     sort_by_bfs(adj, start, N);
        //     fill(visited.begin(), visited.end(), 0);
        //     path = find_longest_path(start, adj, visited);

        //     fill(visited.begin(), visited.end(), 0);
        //     update_visited(adj, path, visited);
        //     // dfs from the start in other direction
        //     from_zero_path = find_longest_path(start, adj, visited);
        //     path = merge_paths_from_same_start(path, from_zero_path);
        //     solutions.push_back(path);

        //     sort(solutions.begin(), solutions.end(), compare_by_size);
        //     if (solutions.size() > 5)
        //     {
        //         solutions.pop_back();
        //     }
        // }
        // longest_path = solutions[0];
        // cout << "Longest path length from normal: " << longest_path.size() << endl;

        // longest_path = {0, 1, 2, 3, 4, 5};
        // for (int node : longest_path)
        // {
        //     visited[node] = 1;
        //     for (Edge e : adj[node])
        //         visited[e.to] += 1;
        // }
        // vector<int> reversed_tail = reverseTail(adj, longest_path, visited);

        // for (int node : reversed_tail)
        // {
        //     cout << node << " ";
        // }
        // cout << endl;

        // for(int i=0;i<5;i++){
        //     //int start = getRandomActiveNode(adj);
        //     start = rand()%N;
        //     sort_by_bfs(adj,start,N);

        //     fill(visited.begin(),visited.end(),0);
        //     path = find_longest_path(start, adj,visited);

        //     fill(visited.begin(),visited.end(),0);
        //     from_zero_path = find_longest_path(start, adj,visited);
        //     path = merge_paths_from_same_start(path,from_zero_path);

        //     if(path.size() > longest_path.size()){
        //         longest_path = path;
        //     }
        // }

        PathValidationResult validation = validate_path_detailed(path, adj);
        cout << validation.error_message << "\n";
    }
    return 0;
}
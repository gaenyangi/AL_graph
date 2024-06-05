#include<iostream>
#include<fstream>
#include<vector>
#include<queue>
#include<limits>
#include<string>

using namespace std;

/*---------------------------------Binary heap class for prim's algorithm---------------------------------*/

class BinaryHeap {
public:
    BinaryHeap(int capacity);
    void insert(int v, int key);
    void decreaseKey(int v, int key);
    int extractMin();
    bool isEmpty() const;
    bool inHeap(int v) const;

private:
    vector<pair<int, int>> heap; // (key, vertex)
    vector<int> pos; // position of vertex in heap
    int size;
    void heapifyDown(int i);
    void heapifyUp(int i);
};

//constructor of Binary heap
BinaryHeap::BinaryHeap(int capacity) : size(0), pos(capacity, -1) {}

//insertion function
void BinaryHeap::insert(int v, int key) {
    heap.push_back({key, v});
    pos[v] = size;
    size++;
    heapifyUp(size - 1);
}

//decrese key function
void BinaryHeap::decreaseKey(int v, int key) {
    int i = pos[v];
    heap[i].first = key;
    heapifyUp(i);
}

//Extract min function
int BinaryHeap::extractMin() {
    if (size == 0) return -1;
    int minVertex = heap[0].second;
    swap(heap[0], heap[size - 1]);
    pos[heap[0].second] = 0;
    pos[minVertex] = -1;
    heap.pop_back();
    size--;
    heapifyDown(0);
    return minVertex;
}

bool BinaryHeap::isEmpty() const {
    return size == 0;
}

bool BinaryHeap::inHeap(int v) const {
    return pos[v] != -1;
}

void BinaryHeap::heapifyDown(int i) {
    int smallest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;
    if (left < size && heap[left].first < heap[smallest].first) {
        smallest = left;
    }
    if (right < size && heap[right].first < heap[smallest].first) {
        smallest = right;
    }
    if (smallest != i) {
        swap(pos[heap[i].second], pos[heap[smallest].second]);
        swap(heap[i], heap[smallest]);
        heapifyDown(smallest);
    }
}

/*
void BinaryHeap::heapifyUp(int i) {
    while (i && heap[i].first < heap[(i - 1) / 2].first) {
        swap(pos[heap[i].second], pos[heap[(i - 1) / 2].second]);
        swap(heap[i], heap[(i - 1) / 2]);
        i = (i - 1) / 2;
    }
}*/

/*-----------------------------MSt class implementing prim's algorithm with heap-------------------------*/

class MST {
public:
    MST(int vertices);
    void addEdge(int u, int v, int w);
    void primMST(int start, const string& output_file);

private:
    int vertices;
    vector<vector<pair<int, int>>> adj;
};

MST::MST(int vertices) : vertices(vertices), adj(vertices) {}

//push edge of node u/v and v/w into adj matrix
void MST::addEdge(int u, int v, int w) {
    adj[u].push_back({v, w});
    adj[v].push_back({u, w});
}

//constructing minimum spanning tree at output.txt with given input
void MST::primMST(int start, const string& output_file) {
    vector<int> key(vertices, numeric_limits<int>::max());
    vector<int> parent(vertices, -1);
    vector<bool> inMST(vertices, false);
    BinaryHeap pq(vertices);

    key[start] = 0;
    pq.insert(start, 0);

    while (!pq.isEmpty()) {
        int u = pq.extractMin();
        inMST[u] = true;

        for (auto& [v, weight] : adj[u]) {
            if (!inMST[v] && weight < key[v]) {
                key[v] = weight;
                parent[v] = u;
                if (pq.inHeap(v)) {
                    pq.decreaseKey(v, weight);
                } else {
                    pq.insert(v, weight);
                }
            }
        }
    }

    ofstream output_file_stream(output_file);
    for (int i = 0; i < vertices; ++i) {
        output_file_stream << i << "\t" << parent[i] << "\n";
    }
    output_file_stream.close();
}

/*--------------------------------- APSP class for Floyd-warshall algorithm--------------------------------*/

class APSP {
public:
    APSP(int vertices);
    void floydWarshall(const vector<vector<int>>& graph, const string& output_file);

private:
    int vertices;
};

APSP::APSP(int vertices) : vertices(vertices) {}

void APSP::floydWarshall(const vector<vector<int>>& graph, const string& output_file) {
    vector<vector<int>> dist = graph;
    vector<vector<int>> pred(vertices, vector<int>(vertices, -1));

    for (int i = 0; i < vertices; ++i) {
        for (int j = 0; j < vertices; ++j) {
            if (i != j && graph[i][j] != numeric_limits<int>::max()) {
                pred[i][j] = i;
            }
        }
    }

    for (int k = 0; k < vertices; ++k) {
        for (int i = 0; i < vertices; ++i) {
            for (int j = 0; j < vertices; ++j) {
                if (dist[i][k] != numeric_limits<int>::max() && dist[k][j] != numeric_limits<int>::max() &&
                    dist[i][k] + dist[k][j] < dist[i][j]) {
                    dist[i][j] = dist[i][k] + dist[k][j];
                    pred[i][j] = pred[k][j];
                }
            }
        }
    }

    ofstream output_file_stream(output_file);

    output_file_stream << "D " << vertices << "\n";
    for (const auto& row : dist) {
        for (int val : row) {
            if (val == numeric_limits<int>::max()) {
                output_file_stream << "INF ";
            } else {
                output_file_stream << val << " ";
            }
        }
        output_file_stream << "\n";
    }

    output_file_stream << "P " << vertices << "\n";
    for (const auto& row : pred) {
        for (int val : row) {
            output_file_stream << (val == -1 ? to_string(-1) : to_string(val+1)) << " ";
        }
        output_file_stream << "\n";
    }

    output_file_stream.close();
}

/*---------------------------------main function---------------------------------*/

int main(int argc, char* argv[]) {

    //check input parameters
    if (argc != 5) {
        cerr << "Usage: " << argv[0] << " <input_file_sp> <input_file_mst> <output_file_sp> <output_file_mst>" << endl;
        return 1;
    }

    string input_file_sp = argv[1];
    string input_file_mst = argv[2];
    string output_file_sp = argv[3];
    string output_file_mst = argv[4];

    // read MST input file
    ifstream mst_input_file(input_file_mst);
    if (!mst_input_file.is_open()) {
        cerr << "Error opening input file: " << input_file_mst << endl;
        return 1;
    }

    // and get values
    int vertices, edges, start_vertex;
    mst_input_file >> vertices >> edges >> start_vertex;

    // construct heap with input undirected graph
    MST mst(vertices);
    for (int i = 0; i < edges; ++i) {
        int u, v, w;
        mst_input_file >> u >> v >> w;
        mst.addEdge(u, v, w);
    }
    mst_input_file.close();

    // call mst class function with input and write the result into output file
    mst.primMST(start_vertex, output_file_mst);

    // read APSP input file
    ifstream sp_input_file(input_file_sp);
    if (!sp_input_file.is_open()) {
        cerr << "Error opening input file: " << input_file_sp << endl;
        return 1;
    }

    int n;
    sp_input_file >> n;

    vector<vector<int>> graph(n, vector<int>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            sp_input_file >> graph[i][j];
        }
    }
    sp_input_file.close();

    // call APSP class function and write the result to output file
    APSP apsp(n);
    apsp.floydWarshall(graph, output_file_sp);

    return 0;
}

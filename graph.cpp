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

// constructor of Binary heap
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
    int minVertex = heap[0].second; // root is min
    swap(heap[0], heap[size - 1]); // swap root with last element
    pos[heap[0].second] = 0;
    pos[minVertex] = -1;
    heap.pop_back();
    size--;
    heapifyDown(0); // and min_heapify again
    return minVertex; // return smallest
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

void BinaryHeap::heapifyUp(int i) {
    while (i && heap[i].first < heap[(i - 1) / 2].first) {
        swap(pos[heap[i].second], pos[heap[(i - 1) / 2].second]);
        swap(heap[i], heap[(i - 1) / 2]);
        i = (i - 1) / 2;
    }
}

/*-----------------------------MSt class implementing prim's algorithm with heap-------------------------*/

class MST {
public:
    MST(int vertices);
    void addEdge(int u, int v, int w);
    void primMST(int start, const string& output_file);

private:
    int vertices;
    vector<vector<pair<int, int>>> adj; //adjacent matrix
};

//constructor
MST::MST(int vertices) : vertices(vertices), adj(vertices) {}

//push edge of node u/v with weight w into adj matrix
void MST::addEdge(int u, int v, int w) {
    adj[u].push_back({v, w});
    adj[v].push_back({u, w}); // add both way cuz it's undirected graph
}

//constructing minimum spanning tree at output.txt with given input
void MST::primMST(int start, const string& output_file) {
    vector<int> key(vertices, 99999); //default MST init with largest value
    vector<int> parent(vertices, -1); // default pred tree
    vector<bool> inMST(vertices, false);
    BinaryHeap pq(vertices);

    key[start] = 0;
    pq.insert(start, 0);

    while (!pq.isEmpty()) {
        int u = pq.extractMin();
        inMST[u] = true;

        for (size_t i = 0; i < adj[u].size(); ++i) {
            int v = adj[u][i].first;
            int weight = adj[u][i].second;

            //check if shorter path exists and swap
            if (!inMST[v] && weight < key[v]) {
                key[v] = weight;
                parent[v] = u;
                //if exist, decrease it's value
                if (pq.inHeap(v)) {
                    pq.decreaseKey(v, weight);
                } else {
                    //if it's new, add it to heap
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

//constructor
APSP::APSP(int vertices) : vertices(vertices) {}

void APSP::floydWarshall(const vector<vector<int>>& graph, const string& output_file) {
    vector<vector<int>> dist = graph; // input adj matrix
    vector<vector<int>> pred(vertices, vector<int>(vertices, -1)); // pred matrix is init as -1

    //init values for pred matrix
    for (int i = 0; i < vertices; ++i) {
        for (int j = 0; j < vertices; ++j) {
            if (i != j && graph[i][j] != 99999) {
                pred[i][j] = i;
            }
        }
    }

    //triply nested for loop of floyd-warshall
    for (int k = 0; k < vertices; ++k) {
        for (int i = 0; i < vertices; ++i) {
            for (int j = 0; j < vertices; ++j) {
                // check whether dist relaxated through k is shorter, and update
                if (dist[i][k] != 99999 && dist[k][j] != 99999 && dist[i][k] + dist[k][j] < dist[i][j]) {
                    dist[i][j] = dist[i][k] + dist[k][j];
                    pred[i][j] = pred[k][j];
                }
            }
        }
    }

    //print the result to output_sp.txt
    ofstream output_file_stream(output_file);

    //shortest path weight matrix D
    output_file_stream << "D " << vertices << "\n";
    for (int i = 0; i < vertices; ++i) {
        for (int j = 0; j < vertices; ++j) {
            int val = dist[i][j];
            if (val == 99999) {
                output_file_stream << "-1 "; //INF
            } else {
                output_file_stream << val << " ";
            }
        }
        output_file_stream << "\n";
    }

    // predecessor matrix P
    output_file_stream << "P " << vertices << "\n";
    for (int i = 0; i < vertices; ++i) {
        for (int j = 0; j < vertices; ++j) {
            int val = pred[i][j];
            output_file_stream << (val == -1 ? to_string(-1) : to_string(val + 1)) << " "; //-1 for NIL
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

    //parse the input arguments
    string input_file_sp = argv[1];
    string input_file_mst = argv[2];
    string output_file_sp = argv[3];
    string output_file_mst = argv[4];

    // open input_mst.txt file
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

    // open input_sp.txt
    ifstream sp_input_file(input_file_sp);
    if (!sp_input_file.is_open()) {
        cerr << "Error opening input file: " << input_file_sp << endl;
        return 1;
    }

    //get the number of lines
    int n;
    sp_input_file >> n;

    // get the data of adj matrix
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

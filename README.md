1. Write a program to sort the elements of an array using Randomized Quick Sort (the
program should report the number of comparisons).
```
#include <iostream>
#include <cstdlib> // For rand()
#include <ctime>   // For seeding rand()

using namespace std;

// Function to swap two elements
void swap(int &a, int &b) {
    int temp = a;
    a = b;
    b = temp;
}

// Partition function with random pivot selection
int randomizedPartition(int arr[], int low, int high, int &comparisons) {
    int randomPivotIndex = low + rand() % (high - low + 1); // Random index
    swap(arr[randomPivotIndex], arr[high]); // Move random pivot to the end
    int pivot = arr[high];
    int i = low - 1;

    for (int j = low; j < high; j++) {
        comparisons++; // Increment comparison counter
        if (arr[j] <= pivot) {
            i++;
            swap(arr[i], arr[j]);
        }
    }
    swap(arr[i + 1], arr[high]);
    return i + 1;
}

// Randomized Quick Sort
void randomizedQuickSort(int arr[], int low, int high, int &comparisons) {
    if (low < high) {
        int pivotIndex = randomizedPartition(arr, low, high, comparisons);
        randomizedQuickSort(arr, low, pivotIndex - 1, comparisons);
        randomizedQuickSort(arr, pivotIndex + 1, high, comparisons);
    }
}

// Driver function
int main() {
    srand(time(0)); // Seed for random number generation

    int n;
    cout << "Enter the number of elements: ";
    cin >> n;

    int arr[n];
    cout << "Enter the elements of the array:\n";
    for (int i = 0; i < n; i++) {
        cin >> arr[i];
    }

    int comparisons = 0;
    randomizedQuickSort(arr, 0, n - 1, comparisons);

    cout << "Sorted array:\n";
    for (int i = 0; i < n; i++) {
        cout << arr[i] << " ";
    }
    cout << "\nTotal number of comparisons: " << comparisons << endl;

    return 0;
}
```


2. Write a program to find the ith smallest element of an array using Randomized Select.

```
#include <iostream>
#include <cstdlib> // For rand() and srand()
#include <ctime>   // For seeding rand()

using namespace std;

// Function to swap two elements
void swap(int &a, int &b) {
    int temp = a;
    a = b;
    b = temp;
}

// Partition function with a random pivot
int randomizedPartition(int arr[], int low, int high) {
    int pivotIndex = low + rand() % (high - low + 1); // Random pivot
    swap(arr[pivotIndex], arr[high]); // Swap random pivot with the last element
    int pivot = arr[high];
    int i = low - 1;

    for (int j = low; j < high; j++) {
        if (arr[j] <= pivot) {
            i++;
            swap(arr[i], arr[j]);
        }
    }
    swap(arr[i + 1], arr[high]);
    return i + 1;
}

// Function to find the ith smallest element using Randomized Select
int randomizedSelect(int arr[], int low, int high, int i) {
    if (low == high) // Base case: only one element
        return arr[low];

    int pivotIndex = randomizedPartition(arr, low, high);

    int k = pivotIndex - low + 1; // Rank of the pivot element

    if (i == k) // Pivot is the ith smallest element
        return arr[pivotIndex];
    else if (i < k) // Look in the left partition
        return randomizedSelect(arr, low, pivotIndex - 1, i);
    else // Look in the right partition
        return randomizedSelect(arr, pivotIndex + 1, high, i - k);
}

int main() {
    srand(time(0)); // Seed the random number generator

    int n, i;
    
    // Take the number of elements in the array as input
    cout << "Enter the number of elements in the array: ";
    cin >> n;

    int arr[n];
    
    // Take the array elements as input
    cout << "Enter the elements of the array: ";
    for (int j = 0; j < n; j++) {
        cin >> arr[j];
    }

    // Take the ith smallest element to find as input
    cout << "Enter the value of i (the ith smallest element to find): ";
    cin >> i;

    // Check if i is valid
    if (i < 1 || i > n) {
        cout << "Invalid value of i. It should be between 1 and " << n << endl;
        return 1;
    }

    // Find the ith smallest element
    int result = randomizedSelect(arr, 0, n - 1, i);
    cout << "The " << i << "th smallest element is: " << result << endl;

    return 0;
}

```



3. Write a program to determine the minimum spanning tree of a graph using Kruskal’s
algorithm.

```
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

// Structure to represent an edge
struct Edge {
    int u, v, weight;
};

// Structure to represent a Disjoint Set (Union-Find)
struct DisjointSet {
    vector<int> parent, rank;

    DisjointSet(int n) {
        parent.resize(n);
        rank.resize(n, 0);
        for (int i = 0; i < n; i++) {
            parent[i] = i;  // Each node is its own parent initially
        }
    }

    // Find function with path compression
    int find(int u) {
        if (parent[u] != u)
            parent[u] = find(parent[u]);
        return parent[u];
    }

    // Union function with union by rank
    void unite(int u, int v) {
        int rootU = find(u);
        int rootV = find(v);

        if (rootU != rootV) {
            if (rank[rootU] > rank[rootV]) {
                parent[rootV] = rootU;
            } else if (rank[rootU] < rank[rootV]) {
                parent[rootU] = rootV;
            } else {
                parent[rootV] = rootU;
                rank[rootU]++;
            }
        }
    }
};

// Function to implement Kruskal's algorithm
void kruskal(int n, vector<Edge>& edges) {
    // Step 1: Sort edges by weight
    sort(edges.begin(), edges.end(), [](Edge& a, Edge& b) {
        return a.weight < b.weight;
    });

    // Step 2: Initialize Disjoint Set
    DisjointSet ds(n);

    int mstWeight = 0;
    vector<Edge> mstEdges;

    // Step 3: Process edges
    for (auto& edge : edges) {
        int u = edge.u, v = edge.v;

        // If including this edge doesn't form a cycle
        if (ds.find(u) != ds.find(v)) {
            ds.unite(u, v);  // Union the sets
            mstEdges.push_back(edge);
            mstWeight += edge.weight;
        }
    }

    // Output the MST
    cout << "Minimum Spanning Tree Edges:\n";
    for (const auto& edge : mstEdges) {
        cout << edge.u << " - " << edge.v << " : " << edge.weight << endl;
    }
    cout << "Total weight of MST: " << mstWeight << endl;
}

int main() {
    int n, e;
    cout << "Enter the number of vertices: ";
    cin >> n;
    cout << "Enter the number of edges: ";
    cin >> e;

    vector<Edge> edges(e);

    cout << "Enter the edges (u, v, weight):\n";
    for (int i = 0; i < e; i++) {
        cin >> edges[i].u >> edges[i].v >> edges[i].weight;
    }

    // Run Kruskal's algorithm to find the MST
    kruskal(n, edges);

    return 0;
}
```


input


```
Enter the number of vertices: 4
Enter the number of edges: 5
Enter the edges (u, v, weight):
0 1 10
0 2 6
0 3 5
1 3 15
2 3 4
```
4. Write a program to implement the Bellman-Ford algorithm to find the shortest paths
from a given source node to all other nodes in a graph.

```
#include <iostream>
#include <vector>
#include <climits> // For INT_MAX

using namespace std;

// Structure to represent an edge
struct Edge {
    int u, v, weight;
};

// Bellman-Ford Algorithm to find shortest paths from source
void bellmanFord(int V, int E, vector<Edge>& edges, int source) {
    // Step 1: Initialize distances from source to all vertices as INFINITY
    vector<int> dist(V, INT_MAX);
    dist[source] = 0;

    // Step 2: Relax all edges V-1 times
    for (int i = 1; i <= V - 1; i++) {
        for (int j = 0; j < E; j++) {
            int u = edges[j].u;
            int v = edges[j].v;
            int weight = edges[j].weight;

            // If distance to v can be shortened by taking edge u-v
            if (dist[u] != INT_MAX && dist[u] + weight < dist[v]) {
                dist[v] = dist[u] + weight;
            }
        }
    }

    // Step 3: Check for negative weight cycles
    for (int i = 0; i < E; i++) {
        int u = edges[i].u;
        int v = edges[i].v;
        int weight = edges[i].weight;

        if (dist[u] != INT_MAX && dist[u] + weight < dist[v]) {
            cout << "Graph contains negative weight cycle\n";
            return;
        }
    }

    // Step 4: Print the distance from the source to each vertex
    cout << "Shortest distances from source " << source << ":\n";
    for (int i = 0; i < V; i++) {
        if (dist[i] == INT_MAX) {
            cout << "Vertex " << i << ": INF\n";
        } else {
            cout << "Vertex " << i << ": " << dist[i] << "\n";
        }
    }
}

int main() {
    int V, E;
    cout << "Enter the number of vertices: ";
    cin >> V;
    cout << "Enter the number of edges: ";
    cin >> E;

    vector<Edge> edges(E);

    cout << "Enter the edges (u, v, weight):\n";
    for (int i = 0; i < E; i++) {
        cin >> edges[i].u >> edges[i].v >> edges[i].weight;
    }

    int source;
    cout << "Enter the source vertex: ";
    cin >> source;

    // Run Bellman-Ford algorithm
    bellmanFord(V, E, edges, source);

    return 0;
}
```

Input

```
Enter the number of vertices: 5
Enter the number of edges: 8
Enter the edges (u, v, weight):
0 1 -1
0 2 4
1 2 3
1 3 2
1 4 2
2 3 5
3 4 -3
Enter the source vertex: 0
```


5. Write a program to implement a B-Tree.

   
```
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

// B-Tree Node structure
struct BTreeNode {
    vector<int> keys;
    vector<BTreeNode*> children;
    bool leaf;
    int t;  // Minimum degree (defines the range for number of keys)

    BTreeNode(int t, bool leaf);
    void insertNonFull(int k);
    void splitChild(int i, BTreeNode* y);
    void traverse();
};

// Constructor to create a new node
BTreeNode::BTreeNode(int t, bool leaf) {
    this->t = t;
    this->leaf = leaf;
}

// B-Tree class
class BTree {
public:
    BTreeNode* root;
    int t;  // Minimum degree

    BTree(int t);
    void insert(int k);
    void traverse();
};

// Constructor to create an empty B-Tree
BTree::BTree(int t) {
    root = nullptr;
    this->t = t;
}

// Insert a key into the B-Tree
void BTree::insert(int k) {
    if (root == nullptr) {
        root = new BTreeNode(t, true);
        root->keys.push_back(k);
    } else {
        if (root->keys.size() == 2 * t - 1) {
            BTreeNode* s = new BTreeNode(t, false);
            s->children.push_back(root);
            s->splitChild(0, root);
            root = s;
        }
        root->insertNonFull(k);
    }
}

// Insert a key into a node that is not full
void BTreeNode::insertNonFull(int k) {
    int i = keys.size() - 1;
    if (leaf) {
        while (i >= 0 && keys[i] > k) {
            i--;
        }
        keys.insert(keys.begin() + i + 1, k);
    } else {
        while (i >= 0 && keys[i] > k) {
            i--;
        }
        i++;
        if (children[i]->keys.size() == 2 * t - 1) {
            splitChild(i, children[i]);
            if (keys[i] < k) {
                i++;
            }
        }
        children[i]->insertNonFull(k);
    }
}

// Split the child of a node
void BTreeNode::splitChild(int i, BTreeNode* y) {
    BTreeNode* z = new BTreeNode(y->t, y->leaf);
    for (int j = 0; j < t - 1; j++) {
        z->keys.push_back(y->keys[j + t]);
    }
    if (!y->leaf) {
        for (int j = 0; j < t; j++) {
            z->children.push_back(y->children[j + t]);
        }
    }
    y->keys.resize(t - 1);
    y->children.resize(t);
    children.insert(children.begin() + i + 1, z);
    keys.insert(keys.begin() + i, y->keys[t - 1]);
}

// Traverse the tree and print keys
void BTreeNode::traverse() {
    int i;
    for (i = 0; i < keys.size(); i++) {
        if (!leaf) {
            children[i]->traverse();
        }
        cout << " " << keys[i];
    }
    if (!leaf) {
        children[i]->traverse();
    }
}

// Traverse the B-Tree
void BTree::traverse() {
    if (root != nullptr) {
        root->traverse();
    }
}

int main() {
    int t;
    cout << "Enter minimum degree of B-Tree: ";
    cin >> t;

    BTree bTree(t);

    while (true) {
        int choice;
        cout << "\nMenu: \n";
        cout << "1. Insert key\n";
        cout << "2. Traverse tree\n";
        cout << "3. Exit\n";
        cout << "Enter your choice: ";
        cin >> choice;

        switch (choice) {
            case 1:
                int key;
                cout << "Enter key to insert: ";
                cin >> key;
                bTree.insert(key);
                break;

            case 2:
                cout << "B-Tree traversal: ";
                bTree.traverse();
                cout << endl;
                break;

            case 3:
                return 0;

            default:
                cout << "Invalid choice\n";
        }
    }

    return 0;
}
```



6. Write a program to implement the Tree Data structure, which supports the following
operations:
a. Insert
b. Search

```
#include <iostream>
using namespace std;

// Node structure for the tree
struct Node {
    int data;
    Node* left;
    Node* right;

    // Constructor to create a new node
    Node(int val) {
        data = val;
        left = nullptr;
        right = nullptr;
    }
};

// Tree class that supports Insert and Search operations
class Tree {
public:
    Node* root;

    Tree() {
        root = nullptr;
    }

    // Function to insert a value into the tree
    void insert(int val) {
        root = insertRec(root, val);
    }

    // Function to search for a value in the tree
    bool search(int val) {
        return searchRec(root, val);
    }

private:
    // Recursive insert function
    Node* insertRec(Node* node, int val) {
        // If the tree is empty, return a new node
        if (node == nullptr) {
            return new Node(val);
        }

        // Otherwise, recur down the tree
        if (val < node->data) {
            node->left = insertRec(node->left, val); // Insert in left subtree
        } else if (val > node->data) {
            node->right = insertRec(node->right, val); // Insert in right subtree
        }

        // Return the unchanged node pointer
        return node;
    }

    // Recursive search function
    bool searchRec(Node* node, int val) {
        // Base case: If the tree is empty or we've found the key
        if (node == nullptr) {
            return false;
        }
        if (node->data == val) {
            return true;
        }

        // Otherwise, search in the left or right subtree
        if (val < node->data) {
            return searchRec(node->left, val); // Search in left subtree
        } else {
            return searchRec(node->right, val); // Search in right subtree
        }
    }
};

int main() {
    Tree tree;
    int n, val;

    // Get the number of values to insert
    cout << "Enter the number of nodes to insert: ";
    cin >> n;

    // Insert elements into the tree
    cout << "Enter " << n << " values to insert into the tree:\n";
    for (int i = 0; i < n; ++i) {
        cin >> val;
        tree.insert(val);
    }

    // Search for a value in the tree
    cout << "\nEnter value to search: ";
    cin >> val;

    if (tree.search(val)) {
        cout << "Value " << val << " found in the tree." << endl;
    } else {
        cout << "Value " << val << " not found in the tree." << endl;
    }

    return 0;
}
```



7. Write a program to search a pattern in a given text using the KMP algorithm.




```

#include <iostream>
#include <vector>
#include <string>
using namespace std;

// Function to compute the LPS (Longest Prefix Suffix) array
void computeLPSArray(const string& pattern, vector<int>& lps) {
    int length = 0; // Length of the previous longest prefix suffix
    lps[0] = 0; // lps[0] is always 0
    int i = 1;

    // Compute the LPS array
    while (i < pattern.length()) {
        if (pattern[i] == pattern[length]) {
            length++;
            lps[i] = length;
            i++;
        } else {
            // Mismatch after length matches
            if (length != 0) {
                length = lps[length - 1];
            } else {
                lps[i] = 0;
                i++;
            }
        }
    }
}

// Function to perform KMP pattern searching
void KMPSearch(const string& text, const string& pattern) {
    int m = pattern.length();
    int n = text.length();

    // Create the LPS array for the pattern
    vector<int> lps(m);
    computeLPSArray(pattern, lps);

    int i = 0; // index for text
    int j = 0; // index for pattern

    while (i < n) {
        if (pattern[j] == text[i]) {
            i++;
            j++;
        }

        if (j == m) {
            cout << "Pattern found at index " << i - j << endl;
            j = lps[j - 1]; // Use the LPS array to avoid unnecessary comparisons
        }
        else if (i < n && pattern[j] != text[i]) {
            // Mismatch after j matches
            if (j != 0) {
                j = lps[j - 1]; // Use the LPS array to avoid unnecessary comparisons
            } else {
                i++;
            }
        }
    }
}

// Main function
int main() {
    string text, pattern;

    // Input text and pattern
    cout << "Enter the text: ";
    getline(cin, text);
    cout << "Enter the pattern: ";
    getline(cin, pattern);

    // Perform KMP search
    KMPSearch(text, pattern);

    return 0;
}
```



8. Write a program to implement a Suffix tree.



```
#include <iostream>
#include <string>
#include <map>
#include <vector>

using namespace std;

// Define the SuffixTreeNode structure
struct SuffixTreeNode {
    // Map for storing child nodes
    map<char, SuffixTreeNode*> children;
    // The index of the suffix that ends at this node
    int startIndex;
    // The length of the suffix (used to store compressed edges)
    int length;

    // Constructor
    SuffixTreeNode() : startIndex(-1), length(0) {}
};

// Define the SuffixTree class
class SuffixTree {
private:
    string text;
    SuffixTreeNode* root;

public:
    // Constructor to initialize the suffix tree with a string
    SuffixTree(const string& str) {
        text = str;
        root = new SuffixTreeNode();
    }

    // Insert a suffix into the suffix tree
    void insertSuffix(const string& suffix, int index) {
        SuffixTreeNode* currentNode = root;

        for (int i = 0; i < suffix.length(); ++i) {
            // Check if the current suffix character already exists in the children
            if (currentNode->children.find(suffix[i]) == currentNode->children.end()) {
                // If not, create a new child node
                SuffixTreeNode* newNode = new SuffixTreeNode();
                newNode->startIndex = index + i;
                newNode->length = suffix.length() - i;
                currentNode->children[suffix[i]] = newNode;
            }

            currentNode = currentNode->children[suffix[i]];
        }
    }

    // Build the suffix tree
    void buildSuffixTree() {
        // Insert all suffixes of the text into the suffix tree
        for (int i = 0; i < text.length(); ++i) {
            insertSuffix(text.substr(i), i);
        }
    }

    // Print the suffix tree (for debugging)
    void printSuffixTree(SuffixTreeNode* node, const string& str) {
        for (auto& child : node->children) {
            string edge = str.substr(child.second->startIndex, child.second->length);
            cout << edge << endl;
            printSuffixTree(child.second, str);
        }
    }

    // Public function to print the tree starting from the root
    void print() {
        printSuffixTree(root, text);
    }
};

int main() {
    string str;
    cout << "Enter a string: ";
    cin >> str;

    // Create the Suffix Tree
    SuffixTree tree(str);

    // Build the tree from the string
    tree.buildSuffixTree();

    // Print the Suffix Tree
    cout << "Suffix Tree:" << endl;
    tree.print();

    return 0;
}
```

//----------------------------------------------------------------------------PREDEFINITIONS--------------------------------------------------------------------------------
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>	// For using the value of pi
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>  // For increasing precision of doubles
#include <limits>   // For setting infinity values
#include <utility>  // For creating pair
#include <queue>    // For using queue

#include <ctime>
#include <algorithm>   
using namespace std;

ofstream out;

//----------------------------------------------------------------------------PREPARATIONS-----------------------------------------------------------------------------------

// Initialize the result values
int result_of_Dijkstra = 0;
int result_of_A_star = 0;
int result_of_A_star_landmark = 0;

// First we need to define the graph G(V,E)
typedef pair <double, double> coordinates; // Define the coordinates of nodes.
typedef pair<int, double> neighbor_and_cost;  // Define the adjacent neighbor node # which is an int and cost which is a double
typedef vector<neighbor_and_cost> set_of_neighbors; // Define set of adjacent neighbors linked to a node

// Then we need to define the type of a vertex
typedef struct new_vertex 
{
	coordinates x_y; // Set coordinates of the node, lat for x, long for y
	set_of_neighbors adjacent_and_number; // Set a vector to store adjacent list of neighbors and costs
} Vertex;

typedef vector<Vertex> graph;	// Record nodes discovered by each algorithm for comparison.

// We need a comparator to make priority queue pop minimum element each time.
struct cost_in_comparision 
{
	bool operator()(const neighbor_and_cost &c1, const neighbor_and_cost &c2) const
	{
		return c1.second > c2.second;
	}
};

// Now we are able to calculate the distance by applying the method in instruction.
double distance(coordinates x_y_1, coordinates x_y_2)
{
	double radius = (2 * M_PI) / 360;   
	double dlat = radius * (x_y_2.first - x_y_1.first);
	double mlat = radius * (x_y_2.first + x_y_1.first) / 2;
	double dlon = radius * (x_y_2.second - x_y_1.second);
	return 6371009 * pow(pow(dlat, 2.0) + pow(cos(mlat)*dlon, 2.0), 0.5);
}

// Input the data
void inputData(graph &G)
{
	int count_of_lines = 1;
	string str = "";
	istringstream input_of_strings;
	ifstream input;
	input.open("graph1000.txt");
	if (!input)
	{
		cerr << " failed to read graph1000.txt" << endl;
		exit(1);
	}

	// We need to access the vertex's number as well as the coordinates
	int iter = 0;
	// Input all of the data
	while (count_of_lines <= G.size())
	{
		input.ignore(10, ':');
		input.ignore(1, ' ');

		getline(input, str, ',');
		double x = atof(str.c_str());
		input.ignore(1);

		getline(input, str);
		double y = atof(str.c_str());

		G[iter].x_y = make_pair(x, y);
		count_of_lines++;
		iter++;
	}

	iter = 0;
	input.ignore(1, '\n');  
	while (iter < 1000)
	{
		input.ignore(10, ':');
		input.ignore(1, ' ');

		getline(input, str);
		input_of_strings.str(str);
		int adjacent_node;
		while (input_of_strings >> adjacent_node)
		{
			if (input_of_strings.peek() == ',')
				input_of_strings.ignore();

			double length = distance(G[iter].x_y, G[adjacent_node - 1].x_y);
			G[iter].adjacent_and_number.push_back(make_pair(adjacent_node - 1, length));
		}

		input_of_strings.clear();

		iter++;
	}
	input.close();  
}

//------------------------------------------------------------------------THREE ALGORITHMS-----------------------------------------------------------------------------------
// 1st:
// Pseudocode of Dijkstra Algorithm with heap
// process: Dijkstra Algorithm(graph G, int start, int end)
// Input: a graph G, an int start indicating the number of the starting node , an int end indicating the number of end node
// Output: the shortestimation path of graph G
/* Step 1:Initialize distances of all vertices as infinite.
   Step 2: Create an empty priority_queue heap. Every item of heap is a pair (distance, vertex). Distance is used used as first item  of pair as first item is by default used to compare two pairs.
   Step 3: Insert source vertex into heap and make its distance as 0.
   Step 4: While either pq doesn't become empty
		a) Extract minimum distance vertex from pq. 
			Let the extracted vertex be u.
		b) Loop through all adjacent of u and do 
			following for every vertex v.
				If there is a shorter path to v
					through u. 
				If dist[v] > dist[u] + weight(u, v)
					Update distance of v, i.e., do
					dist[v] = dist[u] + weight(u, v)
               Insert v into the pq (Even if v is already there)

	thanks to the website: https://www.geeksforgeeks.org/dijkstras-algorithm-for-adjacency-list-representation-greedy-algo-8/
 */
double Dijkstra_Algorithm(const graph &G, const int start, const int end)
{
	int counter = 0;
	ofstream out;
	out.open("First_Part_Dijkstra_Algorithm_Result.txt", ios::app);
	if (!out)
	{
		cerr << "Failed to write to First_Part_Dijkstra_Algorithm_Result.txt" << endl;
		exit(1);
	}

	vector<double> dist(G.size());
	vector<int> pre(G.size());
	vector<int> discovered;
	out << "--------------------------------------------------------------------------------------------------------" << endl;
	out << "Dijkstra Algorithm is in Process:" << endl;
	out << "Starting node = No." << start << endl;
	out << "End node = No." << end << endl;

	for (int i = 0; i < G.size(); i++)
	{
		dist[i] = numeric_limits<double>::infinity(); // First set all distances to infinity
		pre[i] = -1;   // And set no ancestor
	}

	dist[start] = (double)0;

	// Now initialize a priority queue that sets minimum element = highest priority
	priority_queue<neighbor_and_cost, vector<neighbor_and_cost>, cost_in_comparision > Heap;
	Heap.push(make_pair(start, dist[start]));

	while (!Heap.empty())
	{
		int aShort = Heap.top().first; // Find the shortestimation neighbor in the queue
		Heap.pop();

		if (aShort == end)
		{
			out << endl;
			out << "Dijkstra Algorithm's Results: " << endl;
			out << "The shortest distance using Dijkstra Algorithm is " << setprecision(10) << dist[aShort] << " meters." << endl;
			out << "The number of verticies that Dijkstra Algorithm has visited is " << counter << "." << endl;

			result_of_Dijkstra = result_of_Dijkstra + counter;
			out << "--------------------------------------------------------------------------------------------------------" << endl<<endl;
			out.close(); 
			return dist[aShort];
		}

		vector<int>::iterator answer;
		answer = find(discovered.begin(), discovered.end(), aShort);    // To check whether vertex v has been discovered or not
		if (answer == discovered.end()) // When the element has not been discovered yet
		{
			for (int i = 0; i < G[aShort].adjacent_and_number.size(); i++)   // for all edges e of E
			{
				int v = G[aShort].adjacent_and_number[i].first;  // Neighbor's node number
				double c = G[aShort].adjacent_and_number[i].second; // The distance from aShort to its neighbor

				if (dist[v] > dist[aShort] + c)
				{
					dist[v] = dist[aShort] + c;
					pre[v] = aShort;
					Heap.push(make_pair(v, dist[v]));  
				}
			}

			discovered.push_back(aShort);
			counter++;
		}
	
	}
	// When the end node in unreachable
	out << endl;
	out << "Dijkstra Algorithm Underwent Something Wrong. " << endl;
	out << "End node is unreachable!" << endl;
	cout <<"By far, the length is " << dist[end] << endl;
	out << "And the nodes that are discovered are " << counter << endl;
	out << "--------------------------------------------------------------------------------------------------------" << endl << endl;
	out.close(); 
	return dist[end];
}

// 2nd:
// Pseudocode of A* Algorithm
// process: A* Algorithm(graph G, int start,int end)
// Input: a graph G, an int start indicating the number of the starting node , an int end indicating the number of end node
// Output: the shortestimation path of graph G
/*  function reconstruct_path(cameFrom, current)
    total_path := {current}
    while current in cameFrom.Keys:
        current := cameFrom[current]
        total_path.prepend(current)
    return total_path
	// A* finds a path from start to goal.
	// h is the heuristic function. h(n) estimates the cost to reach goal from node n.
	function A_Star(start, goal, h)
    // The set of discovered nodes that may need to be (re-)expanded.
    // Initially, only the start node is known.
    openSet = {start}

    // For node n, cameFrom[n] is the node immediately preceding it on the cheapest path from start to n currently known.
    cameFrom = an empty map

    // For node n, gScore[n] is the cost of the cheapest path from start to n currently known.
    gScore = map with default value of Infinity
    gScore[start] = 0

    // For node n, fScore[n] = gScore[n] + h(n).
    fScore = map with default value of Infinity
    fScore[start] = h(start)

    while openSet is not empty
        current = the node in openSet having the lowest fScore[] value
        if current = goal
            return reconstruct_path(cameFrom, current)

        openSet.Remove(current)
        for each neighbor of current
            // d(current,neighbor) is the weight of the edge from current to neighbor
            // tentative_gScore is the distance from start to the neighbor through current
            tentative_gScore = gScore[current] + d(current, neighbor)
            if tentative_gScore < gScore[neighbor]
                // This path to neighbor is better than any previous one. Record it!
                cameFrom[neighbor] = current
                gScore[neighbor] = tentative_gScore
                fScore[neighbor] = gScore[neighbor] + h(neighbor)
                if neighbor not in openSet
                    openSet.add(neighbor)

    // Open set is empty but goal was never reached
    return failure
	thanks to the website: https://en.wikipedia.org/wiki/A*_search_algorithm#Pseudocode
 */
double A_Star_Algorithm(const graph &G, const int start, const int end)
{
	int counter = 0;
	ofstream out;
	out.open("First_Part_A_Star_Algorithm_Result.txt", ios::app);
	if (!out)
	{
		cerr << "Failed to write to First_Part_A_Star_Algorithm_Result.txt" << endl;
		exit(1);
	}

	vector<double> dist(G.size());
	vector<double> estimation(G.size());
	vector<int> pre(G.size());
	vector<int> discovered;

	out << "--------------------------------------------------------------------------------------------------------" << endl;
	out << "A* Algorithm is in Process:" << endl;
	out << "Starting node = No." << start << endl;
	out << "End node = No." << end << endl;


	for (int i = 0; i < G.size(); i++)
	{
		dist[i] = numeric_limits<double>::infinity(); // Set all distances to infinity
		estimation[i] = numeric_limits<double>::infinity(); // Set all estimationimates to infinity
		pre[i] = -1;   // Set no ancestor
	}

	dist[start] = (double)0;
	double rem_s = distance(G[start].x_y, G[end].x_y);
	estimation[start] = dist[start] + rem_s;
	//out << "The initial estimationimate from start to end is " << setprecision(10) << estimation[start] << endl;

	// Now initialize a priority queue that sets minimum element = highest priority
	priority_queue<neighbor_and_cost, vector<neighbor_and_cost>, cost_in_comparision > Heap;

	Heap.push(make_pair(start, estimation[start]));

	while (!Heap.empty())
	{
		int aShort = Heap.top().first; // Find the shortest neighbor in the queue
		Heap.pop();

		if (aShort == end)
		{
			out << endl;
			out << "A* Algorithm's Results: " << endl;
			out << "The shortest distance using A* Algorithm is " << setprecision(10) << dist[aShort] << " meters." << endl;
			out << "The number of verticies that A* Algorithm has visited is " << counter <<  "." << endl;

			result_of_A_star = result_of_A_star + counter;
			out << "--------------------------------------------------------------------------------------------------------" << endl << endl;
			out.close();
			return dist[aShort];
		}


		vector<int>::iterator answer;
		answer = find(discovered.begin(), discovered.end(), aShort);    // To check whether the vertex v has aleady been discovered or not
		if (answer == discovered.end()) // When the element has not been discovered yet
		{
			for (int i = 0; i < G[aShort].adjacent_and_number.size(); i++)   // for all edges e of E
			{
				int v = G[aShort].adjacent_and_number[i].first;  // Neighbor's node number
				double c = G[aShort].adjacent_and_number[i].second;  // The distance from aShort to its neighbor

				if (dist[v] > dist[aShort] + c)
				{
					dist[v] = dist[aShort] + c;
					pre[v] = aShort;
					double memory_v = distance(G[v].x_y, G[end].x_y);
					estimation[v] = dist[v] + memory_v;
					Heap.push(make_pair(v, estimation[v]));   // Decrease the key
				}
			}
			discovered.push_back(aShort);
			counter++;
		}
	}
	// When the end node in unreachable
	out << endl;
	out << "A* Algorithm Underwent Something Wrong. " << endl;
	out << "End node is unreachable!" << endl;
	cout << "By far, the length is " << dist[end] << endl;
	out << "And the nodes that are discovered are " << counter << endl;
	out << "--------------------------------------------------------------------------------------------------------" << endl << endl;
	out.close();
	return dist[end];
}


// 3rd 
// Pseudocode of A* with landmark Algorithm
// process landmark(graph G, int start, int end, int& random_nodes)
// Input: a graph G, int srart indicating starting node, int end indicating end node, int random_node
// Output: 3 landmarks in random
// Initialize shortest path distances from s to landmarks
// Return the max path of the fore-mentioned shortest path distances

double landmark(const graph &G, int start, int end, const vector<int> &random_nodes)
{
	vector<neighbor_and_cost> memory_landmark;
	vector<neighbor_and_cost> dist_landmark_start;
	vector<neighbor_and_cost> dist_landmark_end;


	// initialize shortestimation path distances from s to landmarks
	for (int i = 0; i < 3; i++)
	{
		int node = random_nodes[i];
		dist_landmark_start.push_back(make_pair(node, distance(G[node].x_y, G[start].x_y)));
		dist_landmark_end.push_back(make_pair(node, distance(G[node].x_y, G[end].x_y)));
		memory_landmark.push_back(make_pair(node, abs(dist_landmark_start[i].second - dist_landmark_end[i].second)));
	}
	return max(max(memory_landmark[0].second, memory_landmark[1].second), memory_landmark[2].second);

}


// process: A* Algorithm with landmark(graph G, int start,int end)
// Input: a graph G, an int start indicating the number of the starting node , an int end indicating the number of end node
// Output: the shortestimation path of graph G
// The rest are similar to the A* Algorithm.
double A_Star_Algorithm_with_Landmarks(const graph &G, const int start, const int end)
{
	int counter = 0;
	ofstream out;
	out.open("First_Part_A_Star_Algorithm_with_Landmarks_Results.txt", ios::app);
	if (!out)
	{
		cerr << "failed to write to First_Part_A_Star_Algorithm_with_Landmarks_Results.txt" << endl;
		exit(1);
	}

	vector<double> dist(G.size());
	vector<double> estimation(G.size());
	vector<int> pre(G.size());
	vector<int> random_nodes(3);
	vector<int> discovered;

	out << "--------------------------------------------------------------------------------------------------------" << endl;
	out << "A* with Landmarks Algorithm is in Process:" << endl;
	out << "Starting node = No." << start << endl;
	out << "End node = No." << end << endl;

	for (int u = 0; u < G.size(); u++)
	{
		dist[u] = numeric_limits<double>::infinity(); // Set all distances to infinity
		estimation[u] = numeric_limits<double>::infinity(); // Set all estimations to infinity
		pre[u] = -1;   // Set no ancestor
	}

	dist[start] = (double)0;

	// Generate three random nodes, and then apply landmark estimations
	for (int i = 0; i < random_nodes.size(); i++)
		random_nodes[i] = rand() % 1000;
	estimation[start] = dist[start] + landmark(G, start, end, random_nodes);

	//out << "The initial estimationimate from start to end is " << setprecision(10) << estimation[start] << endl;

	// Now initialize a priority queue that sets minimum element = highest priority
	priority_queue<neighbor_and_cost, vector<neighbor_and_cost>, cost_in_comparision > Heap;

	Heap.push(make_pair(start, estimation[start]));

	while (!Heap.empty())
	{
		int aShort = Heap.top().first; // the shortestimation neighbor in the queue
		Heap.pop();

		if (aShort == end)
		{
			out << endl;
			out << "A* with Landmarks Algorithm's Results: " << endl;
			out << "The shortest distance using A* with landmarks Algorithm is " << setprecision(10) << dist[aShort] << " meters." << endl;
			out << "The number of verticies that A* with landmarks Algorithm has visited is " << counter <<  "." <<endl;

			result_of_A_star_landmark = result_of_A_star_landmark + counter;
			out << "---------------------------------------------------------------------------------------" << endl << endl;
			out.close(); 						
			return dist[aShort];
		}

		vector<int>::iterator answer;
		answer = find(discovered.begin(), discovered.end(), aShort);    // To check whether the vertex v has aleady been discovered or not
		if (answer == discovered.end()) // When the element has not been discovered yet
		{
			for (int i = 0; i < G[aShort].adjacent_and_number.size(); i++)   // for all edges e of E
			{
				int v = G[aShort].adjacent_and_number[i].first;  // Neighbor node's number
				double c = G[aShort].adjacent_and_number[i].second;  //  The distance from aShort to its neighbor
				if (dist[v] > dist[aShort] + c)
				{
					dist[v] = dist[aShort] + c;
					pre[v] = aShort;
					estimation[v] = dist[v] + landmark(G, v, end, random_nodes);
					Heap.push(make_pair(v, estimation[v]));   // Decrease the key
				}
			}
			discovered.push_back(aShort);
			counter++;
		}
	}
	// When the end node in unreachable
	out << endl;
	out << "A* with Landmarks Algorithm Underwent Something Wrong. " << endl;
	out << "End node is unreachable!" << endl;
	cout << "By far, the length is " << dist[end] << endl;
	out << "And the nodes that are discovered are " << counter << endl;
	out << "--------------------------------------------------------------------------------------------------------" << endl << endl;
	out.close();
	return dist[end];
}

//--------------------------------------------------------------------------MAIN FUNCTION---------------------------------------------------------------------------------
int main()
{
	// Clear
	remove("First_Part_Dijkstra_Algorithm_Result.txt");
	remove("First_Part_A_Star_Algorithm_Result.txt");
	remove("First_Part_A_Star_Algorithm_with_Landmarks_Results.txt");

	// Initialize a Graph containing all 1000 vertices
	int num_of_vertices = 1000;
	graph G(num_of_vertices);
	inputData(G);

	// Initialize random number generator
	srand(time(0));

	// Initialize the distances
	double distance_of_Dijkstra_Algorithm = 0;
	double distance_of_A_Star_Algorithm = 0;
	double distance_of_A_Star_Algorithm_with_Landmarks = 0;

	// Generate 20 random queries for applying Dijkstra, A*, and A*  Landmark
	for (int i = 0; i < 20; i++)
	{
		int random_of_start_node = rand() % 1000;
		int random_of_end_node = rand() % 1000;
		distance_of_Dijkstra_Algorithm += Dijkstra_Algorithm(G, random_of_start_node, random_of_end_node);
		distance_of_A_Star_Algorithm += A_Star_Algorithm(G, random_of_start_node, random_of_end_node);
		distance_of_A_Star_Algorithm_with_Landmarks += A_Star_Algorithm_with_Landmarks(G, random_of_start_node, random_of_end_node);
	}

	double average_node_Dijkstra_Algorithm = ((double)result_of_Dijkstra) / 20.0;
	double average_node_A_Star_Algorithm = ((double)result_of_A_star) / 20.0;
	double average_node_A_Star_Algorithm_with_Landmarks = ((double)result_of_A_star_landmark) / 20.0;
	double average_distance_Dijkstra_Algorithm = distance_of_Dijkstra_Algorithm / 20.0;
	double average_distance_A_Star_Algorithm = distance_of_A_Star_Algorithm / 20.0;
	double average_distance_A_Star_Algorithm_with_Landmarks = distance_of_A_Star_Algorithm_with_Landmarks / 20.0;
	double saving_on_nodes_1 = abs(average_node_A_Star_Algorithm - average_node_Dijkstra_Algorithm);
	double saving_on_nodes_2 = abs(average_node_A_Star_Algorithm_with_Landmarks - average_node_Dijkstra_Algorithm);
	double saving_on_distances_1 = abs(average_distance_A_Star_Algorithm - average_distance_Dijkstra_Algorithm);
	double saving_on_distances_2 = abs(average_distance_A_Star_Algorithm_with_Landmarks - average_distance_Dijkstra_Algorithm);

	ofstream out;
	out.open("Second_Part_Savings.txt");
	out << "--------------------------------------------------------------------------------------------------------" << endl;
	out << "Saving on nodes" << endl;
	out << "Given:" << endl;
	out << "Dijkstra discovered an average of " << average_node_Dijkstra_Algorithm << " nodes." << endl;
	out << "A* discovered an average of " << average_node_A_Star_Algorithm << " nodes." << endl;
	out << "A* with landmarks discovered an average of " << average_node_A_Star_Algorithm_with_Landmarks << " nodes." << endl << endl;
	out << "We can determine that: " << endl;
	out << "When using A* algorithm instead of Dijkstra Algorithm, saving on nodes is " << saving_on_nodes_1 << " nodes." << endl;
	out << "When using A* with landmark algorithm instead of Dijkstra Algorithm, saving on nodes is " << saving_on_nodes_2 << " nodes." << endl<< endl;
	out << "--------------------------------------------------------------------------------------------------------" << endl;
	out << "Saving on distances" << endl;
	out << "Given:" << endl;
	out << "Dijkstra discovered an average of " << average_distance_Dijkstra_Algorithm << " meters." << endl;
	out << "A* discovered an average of " << average_distance_A_Star_Algorithm << " meters." << endl;
	out << "A* with landmarks discovered an average of " << average_distance_A_Star_Algorithm_with_Landmarks << " meters." << endl << endl;
	out << "We can determine that: " << endl;
	out << "When using A* algorithm instead of Dijkstra Algorithm, saving on nodes is " << saving_on_distances_1 << " meters." << endl;
	out << "When using A* with landmark algorithm instead of Dijkstra Algorithm, saving on nodes is " << saving_on_distances_2 << " meters." << endl << endl;
	out << "--------------------------------------------------------------------------------------------------------" << endl;
	out << "We can come to a conclusion that all these algorithms does not change the distance of shortest path, but the total number of vistied nodes." << endl;
	out << "--------------------------------------------------------------------------------------------------------" << endl;
	out << "THIS IS THE END OF THE PROJECT OF CMPT 307-D100: DATA STRUCTURES AND ALGORITHMS." << endl;
	return 0;
}

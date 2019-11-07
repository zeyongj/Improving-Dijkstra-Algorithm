# CMPT-307-Project
This is the introduction of my work of the program project of CMPT 307.

During this project, I first include several header files in the C++ library. 

And then, I make some preparations for reading and handling the data file of graph1000.txt, like defining the definition of graph which is significant in the project for later usage.

After that, I implement the three kinds of algorithms, namely Dijkstra Algorithm, A* Algorithm as well as A* with Landmarks Algorithm to show the shortest path distance and the number of vertices visited each time when using one of these algorithms. The results will be stored in the files with names of “First_Part_XXX_Algorithm.txt”. The pseudocodes of these algorithms are commented in the cpp file, which are taken and inspired from the website: https://www.geeksforgeeks.org/dijkstras-algorithm-for-adjacency-list-representation-greedy-algo-8/ (for Dijkstra Algorithm) as well as the website of https://en.wikipedia.org/wiki/A*_search_algorithm#Pseudocode (for A* Algorithm).

Then, I am able to execute the main function. I set up two variables showing the total nodes and distances for each algorithm, and then call each algorithm for each of 20 queries. Finally, I calculate the average of them by dividing the previous value by 20. As required in the instruction, I also calculate when on average, the savings on nodes and distances when using A* and Landmark Algorithms instead of Dijkstra Algorithm. The result will be stored in the file with name of “Second_Part_Savings.txt”.

According to the result, it is clear that, as for the nodes, when using A* and Landmark Algorithms show significant decrease on the average nodes taken in the process compared to Dijkstra Algorithm. And of course, the shortest distances are all the same.

During the coding work of this project, I got inspired by the website: https://github.com/braeden-mulligan/cmpt307/blob/master/dijkstra.cpp, and some codes are modified from that website. And I would express my great appreciation to the author.

This is the end of the introduction of the project.

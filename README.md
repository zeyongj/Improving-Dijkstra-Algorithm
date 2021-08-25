# CMPT 307 Project

## Introduction
This is the introduction of my work of the program project of CMPT 307.

During this project, I first include several header files in the C++ library. 

And then, I make some preparations for reading and handling the data file of graph1000.txt, like defining the definition of graph which is significant in the project for later usage.

After that, I implement the three kinds of algorithms, namely Dijkstra Algorithm, A* Algorithm as well as A* with Landmarks Algorithm to show the shortest path distance and the number of vertices visited each time when using one of these algorithms. The results will be stored in the files with names of “First_Part_XXX_Algorithm.txt”. The pseudocodes of these algorithms are commented in the cpp file, which are taken and inspired from the website: https://www.geeksforgeeks.org/dijkstras-algorithm-for-adjacency-list-representation-greedy-algo-8/ (for Dijkstra Algorithm) as well as the website of https://en.wikipedia.org/wiki/A*_search_algorithm#Pseudocode (for A* Algorithm).

Then, I am able to execute the main function. I set up two variables showing the total nodes and distances for each algorithm, and then call each algorithm for each of 20 queries. Finally, I calculate the average of them by dividing the previous value by 20. As required in the instruction, I also calculate when on average, the savings on nodes and distances when using A* and Landmark Algorithms instead of Dijkstra Algorithm. The result will be stored in the file with name of “Second_Part_Savings.txt”.

According to the result, it is clear that, as for the nodes, when using A* and Landmark Algorithms show significant decrease on the average nodes taken in the process compared to Dijkstra Algorithm. And of course, the shortest distances are all the same.

During the coding work of this project, I got inspired by the website: https://github.com/braeden-mulligan/cmpt307/blob/master/dijkstra.cpp, and some codes are modified from that website. And I would express my great appreciation to the author.

This is the end of the introduction of the project.

## License

This work is licensed under [Apache License, Version 2.0](https://www.apache.org/licenses/LICENSE-2.0) (or any later version). 

`SPDX-License-Identifier: Apache-2.0-or-later`

## Disclaimer

This repository is ONLY for backup. Students should NEVER use this repository to finish their works, IN ANY WAY.

It is expected that within this course, the highest standards of academic integrity will be maintained, in
keeping with SFU’s Policy S10.01, `Code of Academic Integrity and Good Conduct`.

In this class, collaboration is encouraged for in-class exercises and the team components of the assignments, as well
as task preparation for group discussions. However, individual work should be completed by the person
who submits it. Any work that is independent work of the submitter should be clearly cited to make its
source clear. All referenced work in reports and presentations must be appropriately cited, to include
websites, as well as figures and graphs in presentations. If there are any questions whatsoever, feel free
to contact the course instructor about any possible grey areas.

Some examples of unacceptable behaviour:
- Handing in assignments/exercises that are not 100% your own work (in design, implementation,
wording, etc.), without a clear/visible citation of the source.
- Using another student's work as a template or reference for completing your own work.
- Using any unpermitted resources during an exam.
- Looking at, or attempting to look at, another student's answer during an exam.
- Submitting work that has been submitted before, for any course at any institution.

All instances of academic dishonesty will be dealt with severely and according to SFU policy. This means
that Student Services will be notified, and they will record the dishonesty in the student's file. Students
are strongly encouraged to review SFU’s Code of Academic Integrity and Good Conduct (S10.01) available
online at: http://www.sfu.ca/policies/gazette/student/s10-01.html.

## Author

Zeyong Jin

December 25th, 2019

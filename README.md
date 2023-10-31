# USI AI Cup - TSP Solver
My solution to solve the Traveling Salesman Problem (TSP) as part of the AI course at the USI.
## Table of Contents - cause is cool to have
- [Algorithm](#algorithm)
- [My Results](#my-results)
- [Conclusion](#conclusion)



## Algorithm:
The algorithm used in the code starts with an initial solution generated using the nearest neighbor heuristic. It then employs the 2-opt algorithm to improve this solutionthe simulated annealing as optimization method. The hope is that simulated annealing will allow the algorithm to explore a broader range of solutions and find a tour close to the optimal.

## My Results:
Problem | Best Known | Estimated | Error 
--- | --- | --- | --- 
ch130 |	6110|	6110|	0,00%
d198	|15780|	15852|	0,46%
eil76	|538|	541|	0,56%
fl1577|	22249|	23224|	4,38%
kroa100	|21282	|21343	|0,29%
lin318	|42029|	43600|	3,74%
pcb442	|50778	|52904	|4,19%
pr439	|107217|	112825|	5,23%
rat783|	8806	|9308|	5,70%
u1060	|224094|	241583|	7,80%			
FINAL RESULT| | | 			3%

## Conclusion:
Go check Ant Colony Optimization instead, is better




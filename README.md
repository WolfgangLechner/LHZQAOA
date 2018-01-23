# LHZQAOA

QAOA for fully programmable optimization problem using only nearest neighbor interactions and local fields. The architecture follows Ref. [1] (LHZ).

In QAOA, the unitaries from driver and cost Hamiltonian are:

<a href="https://www.codecogs.com/eqnedit.php?latex=U_x&space;=&space;e^{\beta&space;\sum_{i=1}^N&space;\sigma_x^{(i)}&space;}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?U_x&space;=&space;e^{\beta&space;\sum_{i=1}^N&space;\sigma_x^{(i)}&space;}" title="U_x = e^{\beta \sum_{i=1}^N \sigma_x^{(i)} }" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=U_c&space;=&space;e^{\gamma&space;\sum_{i=1}^N&space;h_i&space;\sigma_z^{(i)}&space;&plus;&space;\sum_{i<j}&space;J_{ij}&space;\sigma_z^{(i)}\sigma_z^{(j)}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?U_c&space;=&space;e^{\gamma&space;\sum_{i=1}^N&space;h_i&space;\sigma_z^{(i)}&space;&plus;&space;\sum_{i<j}&space;J_{ij}&space;\sigma_z^{(i)}\sigma_z^{(j)}}" title="U_c = e^{\gamma \sum_{i<j} J_{ij} \sigma_z^{(i)}\sigma_z^{(j)}}" /></a>

The _challenge_ is that in order to program universal optimization problems the interaction matrix J_ij is all-to-all. A possible solution for that are swap gates to shuffle qubits around which limit the speed of the calculation. As an alternative we propose to make use of the local Hamiltonian Ref. [1] to construct the modified LHZQAOA. The qubit layout is 

<a href="https://www.codecogs.com/eqnedit.php?latex=U_x=&space;\prod_{i=1}^K&space;e^{\beta&space;\sigma_x^{(i)}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?U_x=&space;\prod_{i=1}^K&space;e^{\beta&space;\sigma_x^{(i)}}" title="U_x= \prod_{i=1}^K e^{\beta \sigma_x^{(i)}}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=U_c&space;=&space;\prod_{\{n,o,s,w\}&space;\in&space;p}&space;e^{\Omega&space;\sigma_z^{(n)}&space;\sigma_z^{(o)}\sigma_z^{(s)}\sigma_z^{(w)}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?U_c&space;=&space;\prod_{\{n,o,s,w\}&space;\in&space;p}&space;e^{\Omega&space;\sigma_z^{(n)}&space;\sigma_z^{(o)}\sigma_z^{(s)}\sigma_z^{(w)}}" title="U_c = \prod_{\{n,o,s,w\} \in p} e^{\Omega \sigma_z^{(n)} \sigma_z^{(o)}\sigma_z^{(s)}\sigma_z^{(w)}}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=U_p=&space;\prod_{i=1}^K&space;e^{\gamma&space;\sigma_z^{(i)}&space;}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?U_p=&space;\prod_{i=1}^K&space;e^{\gamma&space;\sigma_z^{(i)}&space;}" title="U_p= \prod_{i=1}^K e^{\gamma \sigma_z^{(i)} }" /></a>


 [1] http://advances.sciencemag.org/content/1/9/e1500838

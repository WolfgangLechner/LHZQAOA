# LHZQAOA

The code demonstrates a fully programmable scheme for QAOA to solve optimization problems using only nearest neighbor interactions and local fields. The qubit architecture follows Ref. [1] (LHZ).

In standard QAOA, the unitaries from driver and cost Hamiltonian are:

<a href="https://www.codecogs.com/eqnedit.php?latex=U_x=&space;\prod_{i=1}^N&space;e^{\beta&space;\sigma_x^{(i)}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?U_x=&space;\prod_{i=1}^N&space;e^{\beta&space;\sigma_x^{(i)}}" title="U_x= \prod_{i=1}^N e^{\beta \sigma_x^{(i)}}" /></a>

and 

<a href="https://www.codecogs.com/eqnedit.php?latex=U_c=&space;\prod_{i<j}&space;e^{\beta&space;J_{ij}\sigma_z^{(i)}\sigma_z^{(j)}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?U_c=&space;\prod_{i<j}&space;e^{\beta&space;J_{ij}\sigma_z^{(i)}\sigma_z^{(j)}}" title="U_c= \prod_{i<j} e^{\beta J_{ij}\sigma_z^{(i)}\sigma_z^{(j)}}" /></a>

The _challenge_ is that in order to program any universal optimization problems the interaction matrix J_ij has to be all-to-all. A possible solution for that are swap gates to shuffle qubits around which limits the speed of the calculation. As an alternative we propose to make use of the local Hamiltonian Ref. [1] to construct the modified LHZQAOA which introduces an quadratic overhead in qubits. The qubit layout is 

![Alt text](img/illustration.png?raw=true "Illustration")

The optimization problem is an all-to-all connected Ising model (left). This is translated via LHZ to a nearest neigbhor model with problem-independent plqauette interactions and local fields that contain the optimization problem (middle). In QAOA, this translates to a model with only nearest neighbor CNOT and RZ gates (right). 

The LHZQAOA contains 3 terms where interactions and programmable local fields can be separated:

<a href="https://www.codecogs.com/eqnedit.php?latex=U_x=&space;\prod_{i=1}^K&space;e^{\beta&space;\sigma_x^{(i)}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?U_x=&space;\prod_{i=1}^K&space;e^{\beta&space;\sigma_x^{(i)}}" title="U_x= \prod_{i=1}^K e^{\beta \sigma_x^{(i)}}" /></a>,

<a href="https://www.codecogs.com/eqnedit.php?latex=U_c&space;=&space;\prod_{\{n,o,s,w\}&space;\in&space;p}&space;e^{\Omega&space;\sigma_z^{(n)}&space;\sigma_z^{(o)}\sigma_z^{(s)}\sigma_z^{(w)}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?U_c&space;=&space;\prod_{\{n,o,s,w\}&space;\in&space;p}&space;e^{\Omega&space;\sigma_z^{(n)}&space;\sigma_z^{(o)}\sigma_z^{(s)}\sigma_z^{(w)}}" title="U_c = \prod_{\{n,o,s,w\} \in p} e^{\Omega \sigma_z^{(n)} \sigma_z^{(o)}\sigma_z^{(s)}\sigma_z^{(w)}}" /></a>,

and 

<a href="https://www.codecogs.com/eqnedit.php?latex=U_p=&space;\prod_{i=1}^K&space;e^{\gamma&space;\sigma_z^{(i)}&space;}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?U_p=&space;\prod_{i=1}^K&space;e^{\gamma&space;\sigma_z^{(i)}&space;}" title="U_p= \prod_{i=1}^K e^{\gamma \sigma_z^{(i)} }" /></a>

The code is a modified version of the MAXCUT example using LHZQAOA with a sequence U = U_x U_c U_p U_x U_c U_p. The FOREST gate-compiler decomposes the constraint unitary U_c with the 4-body interactions to 6 CNOT gates and 1 RZ gate while no(!) ancillas are required. Remarkably, these CNOT gates are also on a flat graph. In the illustration, the red lines indicate the 4-body term in the top plaquette with the sequence: CNOT 1 3, CNOT 3 4, CNOT 4 5, RZ (-2) 5,  CNOT 4 5, CNOT 3 4, CNOT 1 3. The indices in the Figure are the same as in the code. 


 [1] http://advances.sciencemag.org/content/1/9/e1500838

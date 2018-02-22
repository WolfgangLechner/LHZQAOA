##############################################################################
# Copyright 2016-2017 Rigetti Computing
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.
##############################################################################

"""
Finding optimal solution via LHZ with QAOA.
http://advances.sciencemag.org/content/1/9/e1500838
"""

import numpy as np
import pyquil.api as api
from pyquil.paulis import PauliTerm, PauliSum,exponential_map
import networkx as nx
from scipy.optimize import minimize
from grove.pyqaoa.qaoa import QAOA
from lhzqaoa.lhz_qaoa import LHZQAOA


from pyquil.api import get_devices

CXN = api.QVMConnection()
#CXN = api.QPUConnection()

rand_seed = 3
np.random.seed(rand_seed)


#Mapping from n logical qubits to n physical qubits
nmapping = {4: 6, 5: 10, 6: 14, 7: 18}

#For each n logical we use K-N+1 constraints. Labels as on the 19qubit chip
cmappingqpu = {4: [[17, 18, 12],
                  [18, 19, 13],
                  [12, 18, 13, 8]],
              5: [[16, 17, 11],
                  [17, 18, 12],
                  [18, 19, 13],
                  [11, 17, 12, 7],
                  [12, 18, 13, 8],
                  [7, 12, 8, 2]],
              6: [[15, 16, 10],
                  [16, 17, 11],
                  [17, 18, 12],
                  [18, 19, 13],
                  [10, 16, 11, 6],
                  [11, 17, 12, 7],
                  [12, 18, 13, 8],
                  [6, 11, 7, 1],
                  [7, 12, 8, 2]],
              7: [[15, 16, 10],
                  [16, 17, 11],
                  [17, 18, 12],
                  [18, 19, 13],
                  [10, 16, 11, 6],
                  [11, 17, 12, 7],
                  [12, 18, 13, 8],
                  [13, 19, 14, 9],
                  [5, 10, 6, 0],
                  [6, 11, 7, 1],
                  [7, 12, 8, 2]]}

#test for constraint mapping for the qvm with labels from 0 to 5
cmappingqvm = {4: [[0, 1, 3],
                   [1, 2, 4],
                   [1, 3, 4, 5]],
               5: [[0, 1, 4],
                   [1, 2, 5],
                   [2, 3, 6],
                   [1, 4, 5, 7],
                   [2, 5, 6, 8],
                   [5, 7, 8, 9]]
               }

# mapping from jij index to lhz pairs
indextopairmap = {4: [[0, 1], [1, 2], [2, 3], [0, 2], [1, 3], [0, 3]],
                5: [[0, 1], [1, 2], [2, 3], [3, 4], [0, 2], [1, 3], [2, 4], [0, 3], [1, 4], [0, 4]],
                6: [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [0, 2], [1, 3], [2, 4], [3, 5], [0, 3],
                    [1, 4], [2, 5], [0, 4], [1,5], [0, 5]],
                7: [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [0, 2], [1, 3], [2, 4], [3, 5],
                    [4, 6], [0, 3], [1, 4], [2, 5], [3, 6], [0, 4], [1, 5], [2, 6], [0, 5], [1, 6], [0, 6]]
                }



def print_fun(x):
    pass

def LHZ_qaoa(Jij, constraints, steps=1, rand_seed=None, connection=None, samples=None,
                initial_beta=None, initial_gamma=None, initial_omega=None,  minimizer_kwargs=None,
                vqe_option=None):
    """
    LHZ set up method

    :param Jij: the interaction matrix as a list
    :param constraints: list of lists of constraints
    :param steps: (Optional. Default=1) Trotterization order for the QAOA algorithm.
    :param rand_seed: (Optional. Default=None) random seed when beta and gamma angles
        are not provided.
    :param connection: (Optional) connection to the QVM. Default is None.
    :param samples: (Optional. Default=None) VQE option. Number of samples
        (circuit preparation and measurement) to use in operator averaging.
    :param initial_beta: (Optional. Default=None) Initial guess for beta parameters.
    :param initial_gamma: (Optional. Default=None) Initial guess for gamma parameters.
    :param intital_omega: (Optional. Default=None) Initial guess for omega parameters.
    :param minimizer_kwargs: (Optional. Default=None). Minimizer optional arguments.  If None set to
        ``{'method': 'Nelder-Mead', 'options': {'ftol': 1.0e-2, 'xtol': 1.0e-2, 'disp': False}``
    :param vqe_option: (Optional. Default=None). VQE optional arguments.  If None set to
        ``vqe_option = {'disp': print_fun, 'return_all': True, 'samples': samples}``

    """

    constraint_operators = []
    localfield_operator = []
    driver_operators = []

    for c in constraints:
        term = PauliTerm("Z", c[0], 2.0)
        for item in c[1:]:
            term *= PauliTerm("Z", item, 1.0)
        # term += PauliTerm("I", 0, -0.5)
        constraint_operators.append(PauliSum([term]))

    for i in range(len(Jij)):
        localfield_operator.append(PauliSum([PauliTerm("Z", i, Jij[i])]))

    for i in range(len(Jij)):
        driver_operators.append(PauliSum([PauliTerm("X", i, -1.0)]))

    if connection is None:
        connection = CXN

    if minimizer_kwargs is None:
        minimizer_kwargs = {'method': 'Nelder-Mead',
                            'options': {'ftol': 1.0e-2, 'xtol': 1.0e-2,
                                        'disp': False}}
    if vqe_option is None:
        vqe_option = {'disp': print_fun, 'return_all': True,
                      'samples': samples}

    qaoa_inst = LHZQAOA(connection, len(Jij), steps=steps, constraint_ham=constraint_operators, localfield_ham = localfield_operator,
                     ref_hamiltonian=driver_operators, store_basis=True,
                     rand_seed=rand_seed,
                     init_betas=initial_beta,
                     init_gammas=initial_gamma,
                     init_omegas = initial_omega,
                     minimizer=minimize,
                     minimizer_kwargs=minimizer_kwargs,
                     vqe_options=vqe_option)

    return qaoa_inst



def spinglass_qaoa(Jij, indextopairs, nlogic, steps=1, rand_seed=None, connection=None, samples=None,
                initial_beta=None, initial_gamma=None, initial_omega=None,  minimizer_kwargs=None,
                vqe_option=None):
    """
    Spinglass set up

    :param Jij: the interaction matrix as a list
    :param steps: (Optional. Default=1) Trotterization order for the QAOA algorithm.
    :param rand_seed: (Optional. Default=None) random seed when beta and gamma angles
        are not provided.
    :param connection: (Optional) connection to the QVM. Default is None.
    :param samples: (Optional. Default=None) VQE option. Number of samples
        (circuit preparation and measurement) to use in operator averaging.
    :param initial_beta: (Optional. Default=None) Initial guess for beta parameters.
    :param initial_gamma: (Optional. Default=None) Initial guess for gamma parameters.
    :param intital_omega: (Optional. Default=None) Initial guess for omega parameters.
    :param minimizer_kwargs: (Optional. Default=None). Minimizer optional arguments.  If None set to
        ``{'method': 'Nelder-Mead', 'options': {'ftol': 1.0e-2, 'xtol': 1.0e-2, 'disp': False}``
    :param vqe_option: (Optional. Default=None). VQE optional arguments.  If None set to
        ``vqe_option = {'disp': print_fun, 'return_all': True, 'samples': samples}``

    """

    interaction_operators = []
    driver_operators = []

    for i,edge in enumerate(indextopairs):
        term = PauliTerm("Z", edge[0], Jij[i] ) * PauliTerm("Z", edge[1], 1.0)
        interaction_operators.append(PauliSum([term]))

    for i in range(nlogic):
        driver_operators.append(PauliSum([PauliTerm("X", i, -1.0)]))

    if connection is None:
        connection = CXN

    if minimizer_kwargs is None:
        minimizer_kwargs = {'method': 'Nelder-Mead',
                            'options': {'ftol': 1.0e-2, 'xtol': 1.0e-2,
                                        'disp': False}}
    if vqe_option is None:
        vqe_option = {'disp': print_fun, 'return_all': True,
                      'samples': samples}

    qaoa_inst = QAOA(connection, nlogic, steps=steps, cost_ham=interaction_operators,
                     ref_hamiltonian=driver_operators, store_basis=True,
                     rand_seed=rand_seed,
                     init_betas=initial_beta,
                     init_gammas=initial_gamma,
                     minimizer=minimize,
                     minimizer_kwargs=minimizer_kwargs,
                     vqe_options=vqe_option)

    return qaoa_inst


def getmaxprob(states, probs):
    maxprob = 0
    maxstate = 0
    for state, prob in zip(states, probs):
        if prob > maxprob:
            maxprob = prob
            maxstate = state

    return maxstate,maxprob




def getgroundstate(indextopairs,Jij,nlogic):
    mine = 1000.0
    for ci in range(2**nlogic):
        conf = [-1 if int(d) == 0 else 1 for d in str(bin(ci))[2:].zfill(nlogic)]
        e = 0
        for i,edge in enumerate(indextopairs):
            e += Jij[i]*conf[edge[0]]*conf[edge[1]]
            if e < mine:
                mine = e
                groundstate = conf
    return groundstate,mine


if __name__ == "__main__":

    # for device in get_devices():
    #     if device.is_online():
    #         print('Device {} is online'.format(device.name))
    #

    nlogic = 4
    N = nmapping[nlogic]
    constraints = cmappingqvm[nlogic]
    Jij = np.random.random(N)


    print(getgroundstate(indextopairmap[nlogic],Jij,nlogic))


    #QAOA with fully connected spin glass
    #This cannot be implemented but serves as a base comparison
    spinglassinst = spinglass_qaoa(Jij,indextopairmap[nlogic], nlogic, steps=2, rand_seed=rand_seed, samples=None)
    betas, gammas = spinglassinst.get_angles()
    spinglassprobs = spinglassinst.probabilities(np.hstack((betas, gammas)))

    # Comparison with LHZQAOA
    lhzinst = LHZ_qaoa(Jij,constraints, steps=2, rand_seed=rand_seed, samples=None)
    betas, gammas, omegas = lhzinst.get_angles()
    lhzprobs = lhzinst.probabilities(np.hstack((betas, gammas,omegas)))


    lhzmax,lhzmaxprob = getmaxprob(lhzinst.states, lhzprobs)
    sgmax, sgmaxprob = getmaxprob(spinglassinst.states, spinglassprobs)

    print (lhzmax,lhzmaxprob,sgmax,sgmaxprob)


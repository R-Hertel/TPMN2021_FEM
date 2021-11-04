# Calculation of the first derivative of a function using one-dimensional finite elements.
# Exercise TPMN 2021
#
# R. Hertel, Oct. 28th, 2021
# riccardo.hertel@ipcms.unistra.fr

import numpy as np
import matplotlib.pyplot as plt


def nodes(a, b, N):  # create the one-dimensional "mesh"
    nds = np.sort(np.random.uniform(a, b, N-1))
    nds = np.insert(nds, 0, a)
    nds = np.append(nds, b)
    return nds


a = 1; b = 8; N = 50
assert(b > a)
myNodes = nodes(a, b, N)
# print("Node positions:", myNodes)


def element_sizes(nds):  # calculate the size (area) of each element
    return np.diff(nds)


myElementSizes = element_sizes(myNodes)
# print("Finite-element sizes:", myElementSizes)
# print("Sum of element sizes:", sum(myElementSizes))
assert(np.isclose(sum(myElementSizes), b-a))  # np.isclose() checks whether two numbers are nominally equal


def el_nodes(N):  # set up the list of node numbers of each element (two nodes for each element in 1D case)
    nodes_of_elements = []
    for i in range(0, N):
        nodes_of_elements.append((i, i+1))
    return nodes_of_elements


myElementNodes = el_nodes(N)
# print("List of element connectivities:", myElementNodes)


def calc_shape_func_grads(h, N):  # calculate gradient of shape functions (two per element)
    eta_gradients = np.zeros((N, 2))
    for i in range(0, N):
        eta_gradients[i] = (-1/h[i], 1/h[i])
    return eta_gradients


def gradient_in_elements(eta_grads, elementNodes, fval):  # approximate first derivative in element
    func_grad_el = np.zeros(N)
    for i in range(0, N):  # loop over all elements
        for j in range(0, 2):  # loop over nodes of element i
            node = elementNodes[i][j]
            func_grad_el[i] += fval[node] * etaGradients[i, j]
    return func_grad_el


def my_test_function(x):
    f = x * x / 3. + 2.
    return f


def analytic_derivative(x):
    return 2. * x / 3.  # derivative of the test function (for comparison with computed values)


etaGradients = calc_shape_func_grads(myElementSizes, N)
func_vals = []
deriv_at_nodes_exact = []
for i in range(0, N+1):  # define analytic values (function and derivative) at the discretization points
    x = myNodes[i]
    func_vals.append(my_test_function(x))
    deriv_at_nodes_exact.append(analytic_derivative(x))


def calc_node_area(h, elementNodes):  # calculate the area attributed to each node
    node_area = np.zeros(N+1)
    for i in range(0, N):  # loop over all elements
        for j in range(0, 2):  # loop over nodes of element i
            node = elementNodes[i][j]
            node_area[node] += h[i] / 2.
    return node_area


def gradient_at_nodes(func_grad_el, elementNodes, h):
    # calculate weighted average of derivative in neighboring elements
    grad_at_nds = np.zeros(N+1)
    node_weight = calc_node_area(h, elementNodes)
    assert (np.isclose(sum(node_weight), sum(h)))

    for i in range(0, N):  # loop over all elements
        for j in range(0, 2):  # loop over nodes of element i
            node = elementNodes[i][j]
            grad_at_nds[node] += func_grad_el[i] * 0.5 * h[i] / node_weight[node]
    return grad_at_nds


f_gradient_el = gradient_in_elements(etaGradients, myElementNodes, func_vals)
f_grad_nodes = gradient_at_nodes(f_gradient_el, myElementNodes, myElementSizes)


# print("Gradients in elements  = ", f_gradient_el)
# print("Gradients at nodes = ", f_grad_nodes)
# print("Analytic values:", deriv_at_nodes_exact)

# Display results:
plt.plot(myNodes, deriv_at_nodes_exact)
plt.scatter(myNodes, f_grad_nodes, color='red')
plt.show()
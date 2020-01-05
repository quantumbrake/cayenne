Examples
========

Here we discuss some example systems and how to code them up using ``pyssa``.

Zero order system
-----------------
.. math::

    \phi &\xrightarrow[]{k_1} A\\
    \\
    k_1 &= 1.1\\
    A(t=0) &= 100\\

This can be coded up with::

    >>> import numpy as np
    >>> from pyssa import Simulation
    >>> V_r = np.array([[0]])
    >>> V_p = np.array([[1]])
    >>> X0 = np.array([100], dtype=np.int64)
    >>> k_det = np.array([1.1])
    >>> sim = Simulation(V_r, V_p, X0, k_det)
    >>> sim.simulate()
    >>> sim.plot()

.. image:: ../docs/images/ex_0.png
    :scale: 70%
    :align: center
    :alt: Plot of a zero order system.

First order system
------------------

.. math::

    A &\xrightarrow[]{k_1} B\\
    \\
    k_1 &= 1.1\\
    A(t=0) &= 100\\
    B(t=0) &= 20\\

This can be coded up with::

    >>> import numpy as np
    >>> from pyssa import Simulation
    >>> V_r = np.array([[1], [0]])
    >>> V_p = np.array([[0], [1]])
    >>> X0 = np.array([100, 20], dtype=np.int64)
    >>> k_det = np.array([1.1])
    >>> sim = Simulation(V_r, V_p, X0, k_det)
    >>> sim.simulate()
    >>> sim.plot()

.. image:: ../docs/images/ex_1a.png
    :scale: 70%
    :align: center
    :alt: Plot of a first order system.

Suppose you want to use the ``tau_leaping`` algorithm, run 20 repetitions and plot only species :math:`B`. Then do::

    >>> sim.simulate(algorithm="tau_leaping", n_rep=20)
    >>> sim.plot(plot_indices=[1], names=["B"])

.. image:: ../docs/images/ex_1b.png
    :scale: 70%
    :align: center
    :alt: Plot of a first order system with more repetitions.


Enzyme kinetics (second order system with multiple reactions)
-------------------------------------------------------------

.. math::

    \text{Binding}: S + E &\xrightarrow{k1} SE \\
    \text{Dissociation}:SE &\xrightarrow{k2} S + E \\
    \text{Conversion}: SE &\xrightarrow{k3} P + E \\
    \\
    k1 &= 0.006 \\
    k2 &= 0.005 \\
    k3 &= 0.1 \\
    S(t=0) &= 200\\
    E(t=0) &= 50\\
    SE(t=0) &= 0\\
    P(t=0) &= 0\\

This can be coded up with::

    >>> import numpy as np
    >>> from pyssa import Simulation
    >>> V_r = np.array([[1, 0, 0], [1, 0, 0], [0, 1, 1], [0, 0, 0]])
    >>> V_p = np.array([[0, 1, 0], [0, 1, 1], [1, 0, 0], [0, 0, 1]])
    >>> X0 = np.array([200, 50, 0, 0], dtype=np.int64)
    >>> k_det = np.array([0.006, 0.005, 0.1])
    >>> sim = Simulation(V_r, V_p, X0, k_det)
    >>> sim.simulate(max_t=50, n_rep=10)
    >>> sim.plot(names=["S", "E", "SE", "P"])

.. image:: ../docs/images/ex_2a.png
    :scale: 70%
    :align: center
    :alt: Plot of enzyme kinetics simulation.

Since this is a second order system, the size of the system affects the reaction rates. What happens in a larger system? ::

    >>> import numpy as np
    >>> from pyssa import Simulation
    >>> V_r = np.array([[1, 0, 0], [1, 0, 0], [0, 1, 1], [0, 0, 0]])
    >>> V_p = np.array([[0, 1, 0], [0, 1, 1], [1, 0, 0], [0, 0, 1]])
    >>> X0 = np.array([200, 50, 0, 0], dtype=np.int64)
    >>> k_det = np.array([0.006, 0.005, 0.1])
    >>> sim = Simulation(V_r, V_p, X0, k_det, volume=5.0)
    >>> sim.simulate(max_t=50, n_rep=10)
    >>> sim.plot(names=["S", "E", "SE", "P"])

.. image:: ../docs/images/ex_2b.png
    :scale: 70%
    :align: center
    :alt: Plot of enzyme kinetics simulation with a larger volume.

Here we see that the reaction proceeds slower. Less of the product is formed by ``t=50`` compared to the previous case.

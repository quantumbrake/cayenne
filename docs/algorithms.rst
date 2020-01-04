.. _algorithms:

Algorithms
==========

``pyssa`` currently has 3 algorithms:

1. :ref:`direct` (accurate, may be slow)

2. :ref:`tau_leaping` (approximate, faster, needs to be tuned)

3. :ref:`tau_adaptive` (approximate, faster, largely self-tuning)

Methods are described more in depth below.


.. _direct:

Gillespie's direct method (``direct``)
---------------------------------------

.. automodule:: pyssa.algorithms.direct

.. autoclass:: pyssa.algorithms.direct
   :members:
   :undoc-members:
   :show-inheritance:

.. _tau_leaping:

Tau leaping method (``tau_leaping``)
------------------------------------

.. automodule:: pyssa.algorithms.tau_leaping

.. autoclass:: pyssa.algorithms.tau_leaping
   :members:
   :undoc-members:
   :show-inheritance:


.. _tau_adaptive:

Adaptive tau leaping method (experimental, ``tau_adaptive``)
-------------------------------------------------------------

.. automodule:: pyssa.algorithms.tau_adaptive

.. autoclass:: pyssa.algorithms.tau_adaptive
   :members:
   :undoc-members:
   :show-inheritance:

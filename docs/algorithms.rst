.. _algorithms:

Algorithms
==========

``cayenne`` currently has 3 algorithms:

1. :ref:`direct` (accurate, may be slow)

2. :ref:`tau_leaping` (approximate, faster, needs to be tuned)

3. :ref:`tau_adaptive` (approximate, faster, largely self-tuning)

Methods are described more in depth below.


.. _direct:

Gillespie's direct method (``direct``)
---------------------------------------

.. automodule:: cayenne.algorithms.direct

.. autoclass:: cayenne.algorithms.direct
   :members:
   :undoc-members:
   :show-inheritance:

.. _tau_leaping:

Tau leaping method (``tau_leaping``)
------------------------------------

.. automodule:: cayenne.algorithms.tau_leaping

.. autoclass:: cayenne.algorithms.tau_leaping
   :members:
   :undoc-members:
   :show-inheritance:


.. _tau_adaptive:

Adaptive tau leaping method (experimental, ``tau_adaptive``)
-------------------------------------------------------------

.. automodule:: cayenne.algorithms.tau_adaptive

.. autoclass:: cayenne.algorithms.tau_adaptive
   :members:
   :undoc-members:
   :show-inheritance:

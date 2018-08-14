"""This is the docstring for the example.py module.  Modules names should
have short, all-lowercase names.  The module name may have underscores if
this improves readability.

Every module should have a docstring at the very top of the file.  The
module's docstring may extend over multiple lines.  If your docstring does
extend over multiple lines, the closing three quotation marks must be on
a line by itself, preferably preceded by a blank line.

"""
import numpy as np

def direct_naive(V_r, V_p, X0, k_det, max_t = 1, max_iter = 100, volume = 1):
    r"""Naive implementation of the Direct Method.

    Several sentences providing an extended description. Refer to
    variables using back-ticks, e.g. `var`.

    Parameters
    ----------
    var1 : array_like
        Array_like means all those objects -- lists, nested lists, etc. --
        that can be converted to an array.  We can also refer to
        variables like `var1`.
    var2 : int
        The type above can either refer to an actual Python type
        (e.g. ``int``), or describe the type of the variable in more
        detail, e.g. ``(N,) ndarray`` or ``array_like``.
    long_var_name : {'hi', 'ho'}, optional
        Choices in brackets, default first when optional.

    Returns
    -------
    type
        Explanation of anonymous return value of type ``type``.
    describe : type
        Explanation of return value named `describe`.
    out : type
        Explanation of `out`.
    type_without_description

    Other Parameters
    ----------------
    only_seldom_used_keywords : type
        Explanation
    common_parameters_listed_above : type
        Explanation

    Raises
    ------
    BadException
        Because you shouldn't have done that.

    See Also
    --------
    otherfunc : relationship (optional)
    newfunc : Relationship (optional), which could be fairly long, in which
              case the line wraps here.
    thirdfunc, fourthfunc, fifthfunc

    Notes
    -----
    Notes about the implementation algorithm (if needed).

    This can have multiple paragraphs.

    You may include some math:

    .. math:: X(e^{j\omega } ) = x(n)e^{ - j\omega n}

    And even use a Greek symbol like :math:`\omega` inline.

    References
    ----------
    Cite the relevant literature, e.g. [1]_.  You may also cite these
    references in the notes section above.

    .. [1] O. McNoleg, "The integration of GIS, remote sensing,
       expert systems and adaptive co-kriging for environmental habitat
       modelling of the Highland Haggis using object-oriented, fuzzy-logic
       and neural-network techniques," Computers & Geosciences, vol. 22,
       pp. 585-588, 1996.

    Examples
    --------
    These are written in doctest format, and should illustrate how to
    use the function.

    >>> a = [1, 2, 3]
    >>> print [x + 3 for x in a]
    [4, 5, 6]
    >>> print "a\n\nb"
    a
    b

    """

    Na = 6.023e23 # Avogadro's constant
    ite = 1 # Iteration counter
    t = 0 # Time in seconds
    nr = V_r.shape[0] # Number of reactions
    ns = V_r.shape[1] # Number of species
    V = V_p-V_r # nr x ns
    Xt = np.copy(X0) # Number of species at time t
    Xtemp = np.zeros(nr) # Temporary X for updating
    kstoc = np.zeros(nr) # Stochastic rate constants
    orders = np.sum(V_r,1) # Order of rxn = number of reactants

    if np.max(orders) > 3:
        raise RuntimeError('Order greater than 3 detected.')

    if np.max(orders) >1:
        raise RuntimeWarning('Order greater than 1, using volume = ', volume)

    # Determine kstoc from kdet and the highest order or reactions
    for ind in range(nr):
        # If highest order is 3
        if np.max(V_r[ind,:]) == 3:
            kstoc[ind] = k_det[ind] * 6 / np.power(Na*volume, 3)
        elif np.max(V_r[ind,:]) == 2: # Highest order is 2
            kstoc[ind] = k_det[ind] * 2 / np.power(Na*volume, orders[ind])
        else:
            kstoc[ind] = k_det[ind]
    prop = np.copy(kstoc) # Vector of propensities

    while ite<max_iter:
        # Calculate propensities
        for ind1 in range(nr):
            for ind2 in range(ns):
                # prop = kstoc * product of (number raised to order)
                prop[ind1] *= np.power(Xt[ind2], V_r[ind1,ind2])
        # Roulette wheel
        prop0 = np.sum(prop) # Sum of propensities
        prop = prop/prop0 # Normalize propensities to be < 1
        # Concatenate 0 to list of probabilities
        probs = [0]+list(np.cumsum(prop))
        r1 = np.random.rand() # Roll the wheel
        # Identify where it lands and update that reaction
        for ind1 in range(nr+1):
            if r1 <= probs[ind1]:
                Xtemp = Xt + V[ind1-1,:]
                break
        print(Xt, Xtemp, probs)
        print("-------------------------------------------------")
        prop = np.copy(kstoc)
        ite += 1
        # If negative species produced, reject step
        if np.min(Xtemp) < 0:
            continue
        # Update Xt and t
        else:
            Xt = Xtemp
            r2 = np.random.rand()
            t += 1/prop0 * np.log(1/r2)
            if t>max_t:
                print("Reached maximum time (t = )",t)
                break

if __name__ == "__main__":
    V_r = np.array([[1,0,0],[0,1,0]])
    V_p = np.array([[0,1,0],[0,0,1]])
    X0 = np.array([100,0,0])
    k = np.array([1,1])
    print(V_r, X0, V_r.shape, X0.shape)
    direct_naive(V_r, V_p, X0, k, max_t = 10, max_iter = 60)

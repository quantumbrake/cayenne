"""
    Common configuration for all the tests
"""

import csv
import pathlib

import numpy as np
import pytest

from cayenne.utils import Na


CURR_DIR = pathlib.Path(__file__).parent


@pytest.fixture
def setup_basic():
    """
        Setup the basic system.

        Notes
        -----
        A --> B; k = 1.0
        B --> C; k = 1.0

        A0 = 100, B0 = C0 = 0
    """
    V_r = np.array([[1, 0], [0, 1], [0, 0]])
    V_p = np.array([[0, 0], [1, 0], [0, 1]])
    X0 = np.array([100, 0, 0], dtype=np.int64)
    k = np.array([1.0, 1.0])
    species_names = ["A", "B", "C"]
    rxn_names = ["r1", "r2"]
    return species_names, rxn_names, V_r, V_p, X0, k


@pytest.fixture
def setup_large():
    """
        Setup the large system.

        Notes
        -----
        A --> B; k = 1.0
        B --> C; k = 1.0
        C --> D; k = 1.0
        D --> E; k = 1.0
        E --> _; k = 1.0

        A0 = 10, B0 = C0 = D0 = E0 = 0
    """
    V_r = np.array(
        [
            [1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0],
            [0, 0, 1, 0, 0],
            [0, 0, 0, 1, 0],
            [0, 0, 0, 0, 1],
        ]
    )
    V_p = np.array(
        [
            [0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0],
            [0, 0, 1, 0, 0],
            [0, 0, 0, 1, 0],
        ]
    )
    X0 = np.array([10, 0, 0, 0, 0], dtype=np.int64)
    k = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
    species_names = ["A", "B", "C", "D", "E"]
    rxn_names = ["r1", "r2", "r3", "r4", "r5"]
    return species_names, rxn_names, V_r, V_p, X0, k


@pytest.fixture
def setup_system():
    """
        Setup minimal system.

        Contains only ``k_det`` and ``volume``. Used by TestKstoc class as
        filler.
    """
    k_det = np.array([3.0])
    volume = 7.0
    return k_det, volume


@pytest.fixture
def setup_bifurcation():
    """
        Setup bifurcation system.

        Notes
        -----
        A --> B; k = 1.0
        B + D --> B + B; k = 0.01*Na
        A + B --> C; k = 1.0

        A0 = 1, B0 = C0 = 0, D0 = 10
    """
    V_r = np.array([[1, 0, 1], [0, 1, 0], [0, 0, 0], [0, 1, 0]])
    V_p = np.array([[0, 0, 0], [1, 2, 0], [0, 0, 1], [0, 0, 0]])
    k = np.array([1.0, 0.01 * Na, 1.0])
    X0 = np.array([1, 0, 0, 10], dtype=np.int64)
    species_names = ["A", "B", "C", "D"]
    rxn_names = ["r1", "r2", "r3"]
    return species_names, rxn_names, V_r, V_p, X0, k


@pytest.fixture
def setup_long():
    """
        Setup a system that runs for a long time.

        Notes
        -----
        A --> B; k = 1.0
        B --> C; k = 1.0

        A0 = 4e5, B0 = 1000, C0 = 0
    """
    V_r = np.array([[1, 0], [0, 1], [0, 0]])
    V_p = np.array([[0, 0], [1, 0], [0, 1]])
    k = np.array([1.0, 1.0])
    X0 = np.array([int(4e5), 1000, 0], dtype=np.int64)
    species_names = ["A", "B", "C"]
    rxn_names = ["r1", "r2"]
    return species_names, rxn_names, V_r, V_p, X0, k


def read_results(test_id: str):
    """
        Read the simulation results used for accuracy tests.

        Parameters
        ----------
        test_id: str
            The index of the test number to return the results of.

        Returns
        -------
        time_list: List[float]
            Time points at which species amounts are analytically predicted.
        mu_list: List[float]
            Time course of the analytically predicted species amount means.
        std_list: List[float]
            Time course of the analytically predicted species amount standard
            deviations.

        See Also
        --------
        read_results_2sp: Read results for 2 species.

        Notes
        -----
        The accuracy tests are taken from the SBML Test Suite [1]_ .

        References
        ----------
        .. [1] https://github.com/sbmlteam/sbml-test-suite/tree/master/cases/stochastic
    """
    filename = f"{CURR_DIR}/data/results_{test_id}.csv"
    time_list = []
    mu_list = []
    std_list = []
    with open(filename) as fid:
        csv_reader = csv.reader(fid, delimiter=",")
        next(csv_reader)
        for time, mu, std in csv_reader:
            time_list.append(float(time))
            mu_list.append(float(mu))
            std_list.append(float(std))
    return time_list, mu_list, std_list


def read_results_2sp(test_id: str):
    """
        Read the simulation results used for accuracy tests when 2 species
        are tracked.

        Parameters
        ----------
        test_id: str
            The index of the test number to return the results of.

        Returns
        -------
        time_list: List[float]
            Time points at which species amounts are analytically predicted.
        mu_list: List[float]
            Time course of the analytically predicted species amount means.
        std_list: List[float]
            Time course of the analytically predicted species amount standard
            deviations.

        See Also
        --------
        read_results: Read results for one species.

        Notes
        -----
        The accuracy tests are taken from the SBML Test Suite [1]_ .

        References
        ----------
        .. [1] https://github.com/sbmlteam/sbml-test-suite/tree/master/cases/stochastic
    """
    filename = f"{CURR_DIR}/data/results_{test_id}.csv"
    time_list = []
    mu_list = []
    std_list = []
    with open(filename) as fid:
        csv_reader = csv.reader(fid, delimiter=",")
        next(csv_reader)
        for time, mu1, mu2, std1, std2 in csv_reader:
            time_list.append(float(time))
            mu_list.append(np.array([float(mu1), float(mu2)]))
            std_list.append(np.array([float(std1), float(std2)]))
    return time_list, mu_list, std_list


@pytest.fixture
def setup_00001():
    """
        Setup the accuracy test 00001.

        Notes
        -----
        A --> A + A; k = 0.1
        A --> _; k = 0.11

        A0 = 100, max_t = 51, max_iter = 1.5e3
    """
    V_r = np.array([[1, 1]])
    V_p = np.array([[2, 0]])
    X0 = np.array([100], dtype=np.int64)
    k = np.array([0.1, 0.11])
    max_t = 51
    max_iter = int(1.5e3)
    n_rep = 10
    time_list, mu_list, std_list = read_results("00001")
    return (
        V_r,
        V_p,
        X0,
        k,
        time_list,
        np.array(mu_list),
        np.array(std_list),
        max_t,
        max_iter,
        n_rep,
    )


@pytest.fixture
def setup_00003():
    """
        Setup the accuracy test 00003.

        Notes
        -----
        A --> A + A; k = 1.0
        A --> _; k = 1.1

        A0 = 100, max_t = 51, max_iter = 1e5
    """
    species_names = ["A"]
    rxn_names = ["r1", "r2"]
    V_r = np.array([[1, 1]])
    V_p = np.array([[2, 0]])
    X0 = np.array([100], dtype=np.int64)
    k = np.array([1.0, 1.1])
    max_t = 51
    max_iter = int(1e5)
    n_rep = 10
    time_list, mu_list, std_list = read_results("00003")
    return (
        species_names,
        rxn_names,
        V_r,
        V_p,
        X0,
        k,
        time_list,
        np.array(mu_list),
        np.array(std_list),
        max_t,
        max_iter,
        n_rep,
    )


@pytest.fixture
def setup_00004():
    """
        Setup the accuracy test 00004.

        Notes
        -----
        A --> A + A; k = 0.1
        A --> _; k = 0.11

        A0 = 10, max_t = 51, max_iter = 1.5e3
    """
    V_r = np.array([[1, 1]])
    V_p = np.array([[2, 0]])
    X0 = np.array([10], dtype=np.int64)
    k = np.array([0.1, 0.11])
    max_t = 51
    max_iter = int(1.5e3)
    n_rep = 10
    time_list, mu_list, std_list = read_results("00004")
    return (
        V_r,
        V_p,
        X0,
        k,
        time_list,
        np.array(mu_list),
        np.array(std_list),
        max_t,
        max_iter,
        n_rep,
    )


@pytest.fixture
def setup_00005():
    """
        Setup the accuracy test 00005.

        Notes
        -----
        A --> A + A; k = 0.1
        A --> _; k = 0.11

        A0 = 10_000, max_t = 51, max_iter = 5e5
    """
    V_r = np.array([[1, 1]])
    V_p = np.array([[2, 0]])
    X0 = np.array([10000], dtype=np.int64)
    k = np.array([0.1, 0.11])
    max_t = 51
    max_iter = int(5e5)
    n_rep = 10
    time_list, mu_list, std_list = read_results("00005")
    return (
        V_r,
        V_p,
        X0,
        k,
        time_list,
        np.array(mu_list),
        np.array(std_list),
        max_t,
        max_iter,
        n_rep,
    )


@pytest.fixture
def setup_00011():
    """
        Setup the accuracy test 00011.

        Notes
        -----
        A --> A + A; k = 0.1/2
        A --> _; k = 0.11/2

        A0 = 100, max_t = 51, max_iter = 1.5e3
    """
    V_r = np.array([[1, 1]])
    V_p = np.array([[2, 0]])
    X0 = np.array([100], dtype=np.int64)
    # divide k by 2 because rate expression given in units of concentration
    # in the model file
    k = np.array([0.1, 0.11]) / 2
    max_t = 51
    max_iter = int(1.5e3)
    n_rep = 10
    time_list, mu_list, std_list = read_results("00011")
    return (
        V_r,
        V_p,
        X0,
        k,
        time_list,
        np.array(mu_list),
        np.array(std_list),
        max_t,
        max_iter,
        n_rep,
    )


@pytest.fixture
def setup_00020():
    """
        Setup the accuracy test 00020.

        Notes
        -----
        _ --> A; k = 1
        A --> _; k = 0.1

        A0 = 0, max_t = 52, max_iter = 1.5e3
    """
    V_r = np.array([[0, 1]])
    V_p = np.array([[1, 0]])
    X0 = np.array([0], dtype=np.int64)
    k = np.array([1.0, 0.1])
    max_t = 52
    # we did 52 because direct would stop earlier than 50, but we want at least 50 crossed.
    max_iter = int(1.5e3)
    n_rep = 10
    time_list, mu_list, std_list = read_results("00020")
    return (
        V_r,
        V_p,
        X0,
        k,
        time_list,
        np.array(mu_list),
        np.array(std_list),
        max_t,
        max_iter,
        n_rep,
    )


@pytest.fixture
def setup_00021():
    """
        Setup the accuracy test 00021.

        Notes
        -----
        _ --> A; k = 10.0
        A --> _; k = 0.1

        A0 = 0, max_t = 51, max_iter = 1.5e3
    """
    V_r = np.array([[0, 1]])
    V_p = np.array([[1, 0]])
    X0 = np.array([0], dtype=np.int64)
    k = np.array([10.0, 0.1])
    max_t = 51
    max_iter = int(1.5e3)
    n_rep = 10
    time_list, mu_list, std_list = read_results("00021")
    return (
        V_r,
        V_p,
        X0,
        k,
        time_list,
        np.array(mu_list),
        np.array(std_list),
        max_t,
        max_iter,
        n_rep,
    )


@pytest.fixture
def setup_00022():
    """
        Setup the accuracy test 00022.

        Notes
        -----
        _ --> A; k = 5.0
        A --> _; k = 0.1

        A0 = 0, max_t = 51, max_iter = 1.5e3
    """
    V_r = np.array([[0, 1]])
    V_p = np.array([[1, 0]])
    X0 = np.array([0], dtype=np.int64)
    k = np.array([5.0, 0.1])
    max_t = 51
    max_iter = int(1.5e3)
    n_rep = 10
    time_list, mu_list, std_list = read_results("00022")
    return (
        V_r,
        V_p,
        X0,
        k,
        time_list,
        np.array(mu_list),
        np.array(std_list),
        max_t,
        max_iter,
        n_rep,
    )


@pytest.fixture
def setup_00023():
    """
        Setup the accuracy test 00023.

        Notes
        -----
        _ --> A; k = 1000
        A --> _; k = 0.1

        A0 = 0, max_t = 51, max_iter = 1.5e5
    """
    V_r = np.array([[0, 1]])
    V_p = np.array([[1, 0]])
    X0 = np.array([0], dtype=np.int64)
    k = np.array([1000.0, 0.1])
    max_t = 51
    max_iter = int(1.5e5)
    n_rep = 10
    time_list, mu_list, std_list = read_results("00023")
    return (
        V_r,
        V_p,
        X0,
        k,
        time_list,
        np.array(mu_list),
        np.array(std_list),
        max_t,
        max_iter,
        n_rep,
    )


@pytest.fixture
def setup_00030():
    """
        Setup the accuracy test 00030.

        Notes
        -----
        A + A --> A2; k = 0.001
        A2 --> A + A; k = 0.01

        A0 = 100, max_t = 55, max_iter = 1.5e5

        In the model description, they just say k1 = 0.001 without specifying
        deterministic or stochastic. They end up using k1_stoc = 0.001. To have
        k1_stoc = 0.001, we should set k1_det = 0.001/2.
    """
    V_r = np.array([[2, 0], [0, 1]])
    V_p = np.array([[0, 2], [1, 0]])
    X0 = np.array([100, 0], dtype=np.int64)
    k = np.array([0.001 / 2, 0.01])
    max_t = 55
    max_iter = int(1.5e5)
    n_rep = 10
    time_list, mu_list, std_list = read_results_2sp("00030")
    return (
        V_r,
        V_p,
        X0,
        k,
        time_list,
        np.array(mu_list),
        np.array(std_list),
        max_t,
        max_iter,
        n_rep,
    )


@pytest.fixture
def setup_00031():
    """
        Setup the accuracy test 00031.

        Notes
        -----
        A + A --> A2; k = 0.0002
        A2 --> A + A; k = 0.004

        A0 = 1000, max_t = 52, max_iter = 1.5e5

        In the model description, they just say k1 = 0.0002 without specifying
        deterministic or stochastic. They end up using k1_stoc = 0.0002. To
        have k1_stoc = 0.0002, we should set k1_det = 0.0002/2.
    """
    V_r = np.array([[2, 0], [0, 1]])
    V_p = np.array([[0, 2], [1, 0]])
    X0 = np.array([1000, 0], dtype=np.int64)
    k = np.array([0.0002 / 2, 0.004])
    max_t = 52
    max_iter = int(1.5e5)
    n_rep = 10
    time_list, mu_list, std_list = read_results_2sp("00031")
    return (
        V_r,
        V_p,
        X0,
        k,
        time_list,
        np.array(mu_list),
        np.array(std_list),
        max_t,
        max_iter,
        n_rep,
    )


@pytest.fixture
def setup_00037():
    """
        Setup the accuracy test 00037.

        Notes
        -----
        _ --> 5A; k = 1.0
        A --> _; k = 0.2

        A0 = 0, max_t = 53, max_iter = 1.5e3
    """
    V_r = np.array([[0, 1]])
    V_p = np.array([[5, 0]])
    X0 = np.array([0], dtype=np.int64)
    k = np.array([1.0, 0.2])
    max_t = 53
    max_iter = int(1.5e3)
    n_rep = 10
    time_list, mu_list, std_list = read_results("00037")
    return (
        V_r,
        V_p,
        X0,
        k,
        time_list,
        np.array(mu_list),
        np.array(std_list),
        max_t,
        max_iter,
        n_rep,
    )


@pytest.fixture
def setup_00038():
    """
        Setup the accuracy test 00038.

        Notes
        -----
        _ --> 10A; k = 1.0
        A --> _; k = 0.4

        A0 = 0, max_t = 53, max_iter = 1.5e3

        In the model description, they just say k2 = 0.2 without specifying
        deterministic or stochastic. This is a first order reaction, so
        technically k2 = k2_stoc. They instead use kstoc = 0.4 in the model
        file, so this value is used.
    """
    V_r = np.array([[0, 1]])
    V_p = np.array([[10, 0]])
    X0 = np.array([0], dtype=np.int64)
    k = np.array([1.0, 0.4])
    max_t = 53
    max_iter = int(1.5e3)
    n_rep = 10
    time_list, mu_list, std_list = read_results("00038")
    return (
        V_r,
        V_p,
        X0,
        k,
        time_list,
        np.array(mu_list),
        np.array(std_list),
        max_t,
        max_iter,
        n_rep,
    )


@pytest.fixture
def setup_00039():
    """
        Setup the accuracy test 00038.

        Notes
        -----
        _ --> 100A; k = 1.0
        A --> _; k = 4.0

        A0 = 0, max_t = 53, max_iter = 1.5e5

        In the model description, they just say k2 = 0.2 without specifying
        deterministic or stochastic. This is a first order reaction, so
        technically k2 = k2_stoc. They instead use kstoc = 4.0 in the model
        file, so this value is used.
    """
    V_r = np.array([[0, 1]])
    V_p = np.array([[100, 0]])
    X0 = np.array([0], dtype=np.int64)
    k = np.array([1.0, 4.0])
    max_t = 53
    max_iter = int(1.5e5)
    n_rep = 10
    time_list, mu_list, std_list = read_results("00039")
    return (
        V_r,
        V_p,
        X0,
        k,
        time_list,
        np.array(mu_list),
        np.array(std_list),
        max_t,
        max_iter,
        n_rep,
    )


@pytest.fixture
def setup_00001_correct():
    MODEL_00001 = """
        const compartment comp1;
        comp1 = 7; # volume

        r1: S1 => S1 + S1; k1;
        r2: S1 => ; k2;

        k1 = 0.1;
        k2 = 0.11;
        chem_flag = false;

        S1 = 100;
    """
    return MODEL_00001


@pytest.fixture
def setup_00001_nochemflag():
    MODEL_00001 = """
        const compartment comp1;
        comp1 = 7; # volume

        r1: S1 => S1 + S1; k1;
        r2: S1 => ; k2;

        k1 = 0.1;
        k2 = 0.11;

        S1 = 100;
    """
    return MODEL_00001


@pytest.fixture
def setup_00001_norate():
    MODEL_00001 = """
        const compartment comp1;
        comp1 = 7; # volume

        r1: S1 => S1 + S1; k1;
        r2: S1 => ;

        k1 = 0.1;
        chem_flag = false;

        S1 = 100;
    """
    return MODEL_00001


@pytest.fixture
def setup_00001_noratevalue():
    MODEL_00001 = """
        const compartment comp1;
        comp1 = 7; # volume

        r1: S1 => S1 + S1; k1;
        r2: S1 => ; k2;

        k1 = 0.1;
        chem_flag = false;

        S1 = 100;
    """
    return MODEL_00001


@pytest.fixture
def setup_00001_nospeciesvalue():
    MODEL_00001 = """
        const compartment comp1;
        comp1 = 7; # volume

        r1: S1 => S1 + S1; k1;
        r2: S1 => ; k2;

        k1 = 0.1;
        k2 = 0.11;
        chem_flag = false;

    """
    return MODEL_00001


@pytest.fixture
def setup_00001_rateequation():
    MODEL_00001 = """
        const compartment comp1;
        comp1 = 7; # volume

        r1: S1 => S1 + S1; k1*S1;
        r2: S1 => ; k2*S1;

        k1 = 0.1;
        k2 = 0.11;
        chem_flag = false;

        S1 = 100;
    """
    return MODEL_00001


@pytest.fixture
def setup_00001_incompletespeciesvalue():
    MODEL_00001 = """
        const compartment comp1;
        comp1 = 7; # volume

        r1: S1 => S1 + S1; k1;
        r2: S1 => S2; k2;

        k1 = 0.1;
        k2 = 0.11;
        chem_flag = false;

        S1 = 100;

    """
    return MODEL_00001


@pytest.fixture
def setup_00001_nocompartment():
    MODEL_00001 = """
        r1: S1 => S1 + S1; k1;
        r2: S1 => ; k2;

        k1 = 0.1;
        k2 = 0.11;
        chem_flag = false;

        S1 = 100;
    """
    return MODEL_00001

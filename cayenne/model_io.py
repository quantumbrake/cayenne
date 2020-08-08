"""
    The class that handles model IO
"""

import numpy as np

import antimony as sb


class ModelError(Exception):
    pass


class InitialStateError(ModelError):
    pass


class RateConstantError(ModelError):
    pass


class ChemFlagError(ModelError):
    pass


class VolumeError(ModelError):
    pass


class ModelIO:
    """
        Class for loading and parsing models

        Parameters
        ---------
        model_contents : str
            Either the model string or the file path
        content_type : str, {"ModelString", "ModelFile"}
            The type of the model

        Attributes
        ----------
        react_stoic: (ns, nr) ndarray
            A 2D array of the stoichiometric coefficients of the reactants.
            Reactions are columns and species are rows.
        prod_stoic: (ns, nr) ndarray
            A 2D array of the stoichiometric coefficients of the products.
            Reactions are columns and species are rows.
        init_state: (ns,) ndarray
            A 1D array representing the initial state of the system.
        k_det: (nr,) ndarray
            A 1D array representing the deterministic rate constants of the
            system.
        volume: float, optional
            The volume of the reactor vessel which is important for second
            and higher order reactions. Defaults to 1 arbitrary units.
        chem_flag: bool, optional
            If True, divide by Na (Avogadro's constant) while calculating
            stochastic rate constants. Defaults to ``False``.
    """

    def __init__(self, model_contents: str, content_type: str) -> None:
        if content_type == "ModelString":
            er_code = sb.loadAntimonyString(model_contents)
        elif content_type == "ModelFile":
            er_code = sb.loadAntimonyFile(model_contents)
        else:
            raise KeyError(f"Unsupported content_type: {content_type}")
        self.er_code = er_code
        if self.er_code == -1:
            error_msg = "Error while parsing model. Model variable names "
            error_msg += "might be antimony keywords (see docs at https://"
            error_msg += "cayenne.readthedocs.io/en/latest/tutorial.html)."
            raise ModelError(error_msg)
        self.sb_module = sb.getMainModuleName()
        self._parse_model()

    @staticmethod
    def _create_stoic_mat(ns, nr, name_list, stoic_tuple, species_names):
        """ Function to create the stoichiometric matrix """
        stoic_mat = np.zeros([ns, nr], dtype=np.int)
        for index, (names, stoics) in enumerate(zip(name_list, stoic_tuple)):
            for a_name, a_stoic in zip(names, stoics):
                species_index = species_names.index(a_name)
                stoic_mat[species_index, index] += int(a_stoic)
        return stoic_mat

    def _parse_model(self):
        """ Parse model contents """
        react_stoic_tuple = sb.getReactantStoichiometries(self.sb_module)
        react_names = sb.getReactantNames(self.sb_module)
        prod_names = sb.getProductNames(self.sb_module)
        prod_stoic_tuple = sb.getProductStoichiometries(self.sb_module)
        # 0:all, 1:speciescounts, 2:rateconstants, 6:rxnrateequations, 9:compartmentvols
        species_names = sb.getSymbolNamesOfType(self.sb_module, 1)
        self.species_names = species_names
        rxn_names = sb.getSymbolNamesOfType(self.sb_module, 6)
        self.rxn_names = rxn_names
        ns = len(species_names)
        nr = sb.getNumReactions(self.sb_module)

        # Stochastic matrices
        self.react_stoic = self._create_stoic_mat(
            ns, nr, react_names, react_stoic_tuple, species_names
        )
        self.prod_stoic = self._create_stoic_mat(
            ns, nr, prod_names, prod_stoic_tuple, species_names
        )

        # Initial states
        init_state_values = sb.getSymbolInitialAssignmentsOfType(self.sb_module, 1)
        if "" in init_state_values:
            raise InitialStateError("Missing initial value for one of the species.")
        self.init_state = np.array(init_state_values, dtype=np.int64)

        # Rate constants
        rxn_rateeqns = sb.getSymbolEquationsOfType(self.sb_module, 6)
        rxn_rate_names = list(sb.getSymbolNamesOfType(self.sb_module, 2))
        rxn_rate_values = sb.getSymbolInitialAssignmentsOfType(self.sb_module, 2)
        rxn_rate_dict = dict(zip(rxn_rate_names, rxn_rate_values))

        # Chem flag
        try:
            self.chem_flag = True if rxn_rate_dict["chem_flag"] == "true" else False
        except KeyError:
            raise ChemFlagError("The chem flag was not specified in the model.")

        # Check rate constant specifications
        del rxn_rate_dict["chem_flag"]
        rxn_rate_names.remove("chem_flag")
        if "" in rxn_rateeqns:
            raise RateConstantError("Missing rate constant for one of the reactions.")
        for rxn_rateeqn in rxn_rateeqns:
            if rxn_rateeqn not in rxn_rate_names:
                raise RateConstantError(
                    f"{rxn_rateeqn} doesn't match any rate constant."
                )

        # kdet
        try:
            self.k_det = np.array(
                [rxn_rate_dict[rre] for rre in rxn_rateeqns], dtype=float
            )
        except KeyError:
            raise RateConstantError(
                "You are missing a numerical value for one of the rate constants."
            )

        # Volume
        try:
            self.volume = int(
                sb.getSymbolInitialAssignmentsOfType(self.sb_module, 9)[0]
            )
        except IndexError:
            raise VolumeError("Missing compartment information")

    @property
    def args(self):
        """ Returns the attributes of the ModelIO class """
        return (
            self.species_names,
            self.rxn_names,
            self.react_stoic,
            self.prod_stoic,
            self.init_state,
            self.k_det,
            self.chem_flag,
            self.volume,
        )

    @classmethod
    def translate_sbml(cls, sbml_file: str):
        """
            Translate SBML file to Antimony model specification.
            cayenne's model specification is loosely based on Antimony's model
            specification.
        """
        er_code = sb.loadSBMLFile(sbml_file)
        if er_code == -1:
            raise ModelError("Error while parsing model")
        sb_module = sb.getMainModuleName()
        sb_string = sb.getAntimonyString(sb_module)
        return sb_string

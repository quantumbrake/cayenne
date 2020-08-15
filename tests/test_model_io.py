import numpy as np
import pytest

from cayenne.model_io import (
    ChemFlagError,
    InitialStateError,
    ModelError,
    ModelIO,
    RateConstantError,
    VolumeError,
)
from cayenne.simulation import Simulation


class TestModelIO:
    def test_correct_model_str(self, setup_00001_correct, setup_00001):
        modelio = ModelIO(setup_00001_correct, "ModelString")
        (_, _, react_stoic, prod_stoic, init_state, k_det, *_) = modelio.args
        V_r, V_p, X0, k, *_ = setup_00001
        assert (react_stoic == V_r).all()
        assert (prod_stoic == V_p).all()
        assert (X0 == init_state).all()
        assert (k == k_det).all()

    def test_correct_model_file(self, setup_00001):
        modelio = ModelIO("tests/models/00001.txt", "ModelFile")
        (_, _, react_stoic, prod_stoic, init_state, k_det, *_) = modelio.args
        V_r, V_p, X0, k, *_ = setup_00001
        assert (react_stoic == V_r).all()
        assert (prod_stoic == V_p).all()
        assert (X0 == init_state).all()
        assert (k == k_det).all()

    def test_nochemflag(self, setup_00001_nochemflag):
        with pytest.raises(ChemFlagError):
            ModelIO(setup_00001_nochemflag, "ModelString")

    def test_norate(self, setup_00001_norate):
        with pytest.raises(RateConstantError):
            ModelIO(setup_00001_norate, "ModelString")

    def test_noratevalue(self, setup_00001_noratevalue):
        with pytest.raises(RateConstantError):
            ModelIO(setup_00001_noratevalue, "ModelString")

    def test_rateequation(self, setup_00001_rateequation):
        with pytest.raises(RateConstantError):
            ModelIO(setup_00001_rateequation, "ModelString")

    def test_nospeciesvalue(self, setup_00001_nospeciesvalue):
        with pytest.raises(InitialStateError):
            ModelIO(setup_00001_nospeciesvalue, "ModelString")

    def test_incompletespeciesvalue(self, setup_00001_incompletespeciesvalue):
        with pytest.raises(InitialStateError):
            ModelIO(setup_00001_incompletespeciesvalue, "ModelString")

    def test_nocompartment(self, setup_00001_nocompartment):
        with pytest.raises(VolumeError):
            ModelIO(setup_00001_nocompartment, "ModelString")

    def test_translate_sbml(self):
        sb_string = ModelIO.translate_sbml("tests/models/00001.xml")
        assert isinstance(sb_string, str)

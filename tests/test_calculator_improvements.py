from __future__ import annotations

import pytest
from simpyson.calculator import SimpCalc


def test_validation_spinsys():
    """Test that spinsys cannot be None."""
    with pytest.raises(ValueError):
        SimpCalc(spinsys=None)


def test_validation_pulse_sequence():
    """Test invalid pulse sequence types."""
    with pytest.raises(ValueError):
        SimpCalc(spinsys="spinsys { channels 1H }", pulse_sequence=123)


def test_missing_parameters():
    """Test missing required parameters."""
    calc = SimpCalc(spinsys="spinsys { channels 1H }", pulse_sequence="pulse_90")
    with pytest.raises(ValueError, match="Missing required parameters"):
        calc.generate_par()


def test_dry_run():
    """Test dry_run option."""
    calc = SimpCalc(
        spinsys="spinsys { channels 1H }",
        pulse_sequence="pulse_90",
        proton_frequency=400e6,
        spin_rate=10000,
        start_operator="I1z",
        detect_operator="I1p",
        np=1024,
        sw=20000,
        method="direct",
        crystal_file="rep100",
        gamma_angles=10,
        verbose=0,
    )
    # Should not raise FileNotFoundError even if simpson is missing
    cmd = calc.run(dry_run=True)
    assert isinstance(cmd, str)
    assert "simpson" in cmd

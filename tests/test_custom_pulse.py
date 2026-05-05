from __future__ import annotations

import re

from simpyson.calculator import SimpCalc


def test_custom_pulse_sequence_params():
    """Test that custom pulse sequence receives parameters from SimpCalc."""
    code = """
    pulse $par(my_param) 0 0 0 0
    """

    calc = SimpCalc(
        spinsys="spinsys { channels 1H }",
        pulse_sequence=code,
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
        my_param=10.0,       # This must pass to CustomPulseSequence
        unused_param=20.0,   # This should not
    )

    assert 'variable_my_param' in calc.pulse_sequence.parameters
    assert calc.pulse_sequence.parameters['variable_my_param'] == 10.0
    assert 'variable_unused_param' not in calc.pulse_sequence.parameters

    # Verify generation
    par_block = calc.generate_par()
    assert re.search(r"variable\s+my_param\s+10\.0", par_block)


def test_standard_params_in_custom_code():
    """Test using standard parameters in custom code."""
    code = """
    delay $par(np)
    """

    calc = SimpCalc(
        spinsys="spinsys { channels 1H }",
        pulse_sequence=code,
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

    # np is required by code
    assert 'variable_np' in calc.pulse_sequence.parameters

    # It should appear in par block as a standard parameter
    par_block = calc.generate_par()
    assert "np                   1024" in par_block

    # It should not appear as "variable np"
    assert "variable np" not in par_block

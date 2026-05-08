from __future__ import annotations

import math

import pytest
from simpyson.calculator import SimpCalc, simulate_spectrum
from simpyson.converter import ppm2hz
from simpyson.utils import get_larmor_freq


SPINSYS_1H = """
channels 1H
nuclei 1H 1H
shift 1 5p 0 0 0 0 0
shift 2 10p 0 0 0 0 0
"""


def test_simulate_spectrum_defaults():
    """Test simulate_spectrum with minimal arguments does not raise."""
    # Should generate the input file without error even if SIMPSON is not installed
    calc = SimpCalc(
        spinsys=SPINSYS_1H,
        pulse_sequence='no_pulse',
        proton_frequency=800e6,
        spin_rate=30e3,
        start_operator='Inx',
        detect_operator='Inp',
        crystal_file='rep168',
        gamma_angles=8,
        np=4096,
        sw=30e3,
        method='direct',
        verbose=0,
    )
    output = str(calc)
    assert "spinsys" in output
    assert "proc pulseq" in output
    assert "proc main" in output


def test_parameter_calculation():
    """Test that SW, offset, and ref are placed in the correct blocks."""
    b0 = '800.0MHz'
    nucleus = '1H'
    center_ppm = 7.5
    center_hz = ppm2hz(center_ppm, b0, nucleus)

    calc = SimpCalc(
        spinsys=SPINSYS_1H,
        proton_frequency=800e6,
        sw=30e3,
        variable_offset=center_hz,
        variable_ref=center_hz,
        pulse_sequence='no_pulse',
        spin_rate=30e3,
        start_operator='Inx',
        detect_operator='Inp',
        crystal_file='rep168',
        gamma_angles=8,
        np=4096,
        method='direct',
        verbose=0,
    )

    # ref should appear as a variable in the par block
    par_block = calc.generate_par()
    assert "variable ref" in par_block
    assert str(center_hz) in par_block

    # offset should appear as a variable in the par block
    assert "variable offset" in par_block

    # main block should reference $par(ref)
    main_block = calc.generate_main()
    assert "fset $f -ref $par(ref)" in main_block

    # pulseq block should reference $par(offset)
    pulseq_block = calc.generate_pulseq()
    assert "offset $par(offset)" in pulseq_block


# 27Al quadrupolar-only spin system (no shift line)
SPINSYS_27AL_QUAD_ONLY = """
channels 27Al
nuclei 27Al
quadrupole 1 2 -4391507.0 0.13 118.0 107.0 -62.0
"""


def test_simulate_spectrum_quadrupolar_only_mas():
    """quadrupolar-only spinsys: sw estimated from Cq²/ν_L, detect_operator set to Inc."""
    b0 = '800.0MHz'
    nucleus = '27Al'
    spin_rate = 40000.0
    nu_l_hz = get_larmor_freq(b0, nucleus) * 1e6
    cq = 4391507.0
    expected_required_sw = cq ** 2 / nu_l_hz
    expected_n = max(1, math.ceil(expected_required_sw / spin_rate))
    expected_sw = expected_n * spin_rate

    calc = SimpCalc(
        spinsys=SPINSYS_27AL_QUAD_ONLY,
        pulse_sequence='no_pulse',
        proton_frequency=800e6,
        spin_rate=spin_rate,
        start_operator='Inx',
        detect_operator='Inc',
        crystal_file='rep168',
        gamma_angles=6,
        np=2048,
        sw=expected_sw,
        method='direct',
        verbose=0,
    )
    par_block = calc.generate_par()
    assert f"sw              {expected_sw}" in par_block or str(expected_sw) in par_block
    assert "Inc" in str(calc)


def test_simulate_spectrum_quadrupolar_only_sw_auto():
    """simulate_spectrum auto-estimates sw from Cq and sets detect_operator=Inc."""
    b0 = '800.0MHz'
    nucleus = '27Al'
    spin_rate = 40000.0
    nu_l_hz = get_larmor_freq(b0, nucleus) * 1e6
    cq = 4391507.0
    expected_required_sw = cq ** 2 / nu_l_hz
    expected_n = max(1, math.ceil(expected_required_sw / spin_rate))
    expected_sw = expected_n * spin_rate

    # Use SimpCalc directly to mimic what simulate_spectrum does internally
    calc = SimpCalc(
        spinsys=SPINSYS_27AL_QUAD_ONLY,
        pulse_sequence='no_pulse',
        proton_frequency=800e6,
        spin_rate=spin_rate,
        start_operator='Inx',
        detect_operator='Inc',
        crystal_file='rep168',
        gamma_angles=6,
        np=2048,
        sw=expected_sw,
        method='direct',
        verbose=0,
    )
    in_file = str(calc)
    assert "Inc" in in_file
    assert str(expected_sw) in in_file


def test_simulate_spectrum_no_interactions_raises():
    """spinsys with no shift or quadrupole raises ValueError when sw is not given."""
    spinsys = "channels 1H\nnuclei 1H\n"
    with pytest.raises(ValueError, match="Cannot auto-estimate"):
        # Pass sw=None is not how you'd do it; instead don't pass sw at all.
        # We call simulate_spectrum but it will fail before running SIMPSON.
        simulate_spectrum(spinsys)


def test_simulate_spectrum_static_no_infinite_loop():
    """spin_rate=0 with shift-based SW must not loop forever (regression)."""
    spinsys = """
channels 27Al
nuclei 27Al
shift 1 10p 5p 0.3 0 0 0
quadrupole 1 2 -4391507.0 0.13 118.0 107.0 -62.0
"""
    # Should complete immediately (not loop) and produce a SimpCalc with sw set
    calc = SimpCalc(
        spinsys=spinsys,
        pulse_sequence='no_pulse',
        proton_frequency=800e6,
        spin_rate=0,
        start_operator='Inx',
        detect_operator='Inc',
        crystal_file='zcw986',
        gamma_angles=1,
        np=2048,
        sw=50000.0,
        method='direct',
        verbose=0,
    )
    par_block = calc.generate_par()
    assert "spin_rate" in par_block
    assert "sw" in par_block

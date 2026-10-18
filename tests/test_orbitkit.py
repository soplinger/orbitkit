import math
from pathlib import Path

import pytest

import orbitkit
from orbitkit import (
    BODIES,
    bielliptic,
    circular_speed,
    escape_speed,
    fmt_duration,
    hohmann,
    kepler_solve,
    main,
    orbital_period,
    orbital_state,
    plan_mission,
    plane_change_dv,
    rocket_propellant,
)

EXAMPLE = Path(__file__).resolve().parent.parent / "examples" / "leo_to_geo.orbit"
MU_EARTH = BODIES["earth"]["mu"]
R_EARTH = BODIES["earth"]["radius"]


def test_circular_speed_at_200_km_is_about_7_79_km_s():
    r = R_EARTH + 200
    assert circular_speed(MU_EARTH, r) == pytest.approx(7.785, abs=0.01)


def test_period_at_200_km_is_about_88_5_minutes():
    period_min = orbital_period(MU_EARTH, R_EARTH + 200) / 60
    assert 88.0 < period_min < 89.0


def test_escape_speed_at_surface_is_about_11_18_km_s():
    assert escape_speed(MU_EARTH, R_EARTH) == pytest.approx(11.18, abs=0.02)


def test_escape_is_root_two_times_circular():
    r = R_EARTH + 500
    assert escape_speed(MU_EARTH, r) == pytest.approx(math.sqrt(2) * circular_speed(MU_EARTH, r))


def test_leo_to_geo_hohmann_total_is_about_3_93_km_s():
    result = hohmann(MU_EARTH, R_EARTH + 200, R_EARTH + 35786)
    assert 3.92 < result["total_kms"] < 3.94
    assert result["dv1_kms"] > 0 and result["dv2_kms"] > 0


def test_hohmann_transfer_time_is_half_the_transfer_ellipse_period():
    r1, r2 = R_EARTH + 200, R_EARTH + 35786
    result = hohmann(MU_EARTH, r1, r2)
    a_t = (r1 + r2) / 2
    assert result["time_s"] == pytest.approx(orbital_period(MU_EARTH, a_t) / 2)


def test_lowering_orbit_uses_two_retrograde_burns():
    result = hohmann(1.0, 7000.0, 6700.0)
    assert result["dv1_kms"] < 0
    assert result["dv2_kms"] < 0
    assert result["total_kms"] == pytest.approx(abs(result["dv1_kms"]) + abs(result["dv2_kms"]))


def test_hohmann_same_orbit_is_rejected():
    with pytest.raises(ValueError, match="same"):
        hohmann(MU_EARTH, 7000.0, 7000.0)


def test_bielliptic_beats_hohmann_for_large_ratio():
    # Normalised units (mu = 1): radius ratio 15 with apex at 50.
    hoh = hohmann(1.0, 1.0, 15.0)
    bi = bielliptic(1.0, 1.0, 15.0, 50.0)
    assert bi["total_kms"] < hoh["total_kms"]
    assert bi["total_kms"] == pytest.approx(0.5304, abs=1e-3)


def test_bielliptic_apex_must_be_above_both_orbits():
    with pytest.raises(ValueError, match="above both"):
        bielliptic(MU_EARTH, 7000.0, 42000.0, 30000.0)


def test_plane_change_at_90_degrees():
    assert plane_change_dv(7.5, 90) == pytest.approx(2 * 7.5 * math.sin(math.radians(45)))


@pytest.mark.parametrize("deg", [0, -5, 181])
def test_plane_change_rejects_out_of_range_angles(deg):
    with pytest.raises(ValueError):
        plane_change_dv(7.0, deg)


def test_kepler_solution_satisfies_kepler_equation():
    e, m = 0.5, 1.0
    big_e = kepler_solve(m, e)
    assert big_e - e * math.sin(big_e) - m == pytest.approx(0.0, abs=1e-9)


def test_kepler_circular_orbit_has_E_equal_M():
    assert kepler_solve(1.2, 0.0) == pytest.approx(1.2)


def test_state_at_perigee_is_a_times_one_minus_e():
    s = orbital_state(MU_EARTH, 10000.0, 0.5, 0.0)
    assert s["r_km"] == pytest.approx(5000.0)
    assert s["nu_deg"] == pytest.approx(0.0, abs=1e-9)


def test_state_on_circular_orbit_matches_circular_speed():
    s = orbital_state(MU_EARTH, 7000.0, 0.0, 45.0)
    assert s["r_km"] == pytest.approx(7000.0)
    assert s["nu_deg"] == pytest.approx(45.0)
    assert s["v_kms"] == pytest.approx(circular_speed(MU_EARTH, 7000.0))


def test_rocket_propellant_zero_dv_needs_no_propellant():
    assert rocket_propellant(0.0, 300, 1000) == 0


def test_rocket_propellant_for_mass_ratio_two_equals_dry_mass():
    dv_kms = 9.80665 * 300 * math.log(2) / 1000
    assert rocket_propellant(dv_kms, 300, 1000) == pytest.approx(1000.0)


def test_rocket_rejects_negative_dv():
    with pytest.raises(ValueError):
        rocket_propellant(-1.0, 300, 1000)


def test_example_mission_totals():
    report = plan_mission(EXAMPLE.read_text(encoding="utf-8"))
    assert report["errors"] == []
    assert [leg["line"] for leg in report["legs"]] == [4, 5, 6]
    assert 3.93 < report["total_kms"] < 3.94
    assert report["elapsed_s"] == pytest.approx(18933 + 7200, abs=60)
    assert 3700 < report["propellant_kg"] < 3800
    assert report["within_budget"] is True


def test_mission_over_budget_is_flagged():
    text = "body earth\norbit circular 200\nhohmann 35786\nbudget 1\n"
    report = plan_mission(text)
    assert report["within_budget"] is False


def test_errors_carry_line_numbers():
    text = "body earth\norbit circular 200\nhohmann abc\n"
    report = plan_mission(text)
    assert report["errors"] == [{"line": 3, "message": "altitude 'abc' is not a number"}]
    assert report["total_kms"] is None


def test_unknown_directive_points_at_its_line():
    report = plan_mission("body earth\norbit circular 200\nteleport 5\n")
    assert report["errors"][0]["line"] == 3
    assert "unknown directive" in report["errors"][0]["message"]


def test_manoeuvre_before_orbit_is_an_error():
    report = plan_mission("body earth\nhohmann 35786\n")
    assert "orbit" in report["errors"][0]["message"]


def test_body_after_manoeuvre_is_an_error():
    report = plan_mission("body earth\norbit circular 200\nburn 0.1\nbody mars\n")
    assert report["errors"][0]["line"] == 4


def test_missing_body_is_reported_on_line_zero():
    report = plan_mission("burn 0.1\n")
    assert report["errors"] == [{"line": 0, "message": "mission has no 'body' directive"}]


def test_comments_and_blank_lines_are_ignored():
    text = "# plan\n\nbody moon   # the moon\norbit circular 100\n\nburn 0.05 # small\n"
    report = plan_mission(text)
    assert report["errors"] == []
    assert report["legs"][0]["line"] == 6
    assert report["total_kms"] == pytest.approx(0.05)


def test_negative_burn_is_rejected_with_line():
    report = plan_mission("body earth\norbit circular 200\nburn -1\n")
    assert report["errors"][0]["line"] == 3


def test_fmt_duration_formats_hours_and_minutes():
    assert fmt_duration(3725) == "1 h 02 min"


def test_unknown_body_is_rejected():
    with pytest.raises(ValueError, match="unknown body"):
        orbitkit.get_body("pluto")


def test_cli_info_prints_speed(capsys):
    assert main(["info", "earth", "200"]) == 0
    assert "speed:" in capsys.readouterr().out


def test_cli_transfer_with_via_reports_cheaper_option(capsys):
    assert main(["transfer", "earth", "200", "35786", "--via", "100000"]) == 0
    out = capsys.readouterr().out
    assert "cheaper:" in out


def test_cli_plan_example_exits_zero(capsys):
    assert main(["plan", str(EXAMPLE)]) == 0
    assert "Total delta-v" in capsys.readouterr().out


def test_cli_plan_budget_breach_exits_two(tmp_path, capsys):
    mission = tmp_path / "tight.orbit"
    mission.write_text("body earth\norbit circular 200\nhohmann 35786\nbudget 1\n", encoding="utf-8")
    assert main(["plan", str(mission)]) == 2
    assert "OVER BUDGET" in capsys.readouterr().out


def test_cli_plan_error_exits_one(tmp_path, capsys):
    mission = tmp_path / "bad.orbit"
    mission.write_text("body earth\norbit circular 200\nhohmann abc\n", encoding="utf-8")
    assert main(["plan", str(mission)]) == 1
    assert "3: error:" in capsys.readouterr().out


def test_cli_json_output_is_parseable(capsys):
    assert main(["plan", str(EXAMPLE), "--json"]) == 0
    import json

    data = json.loads(capsys.readouterr().out)
    assert data["body"] == "earth"
    assert len(data["legs"]) == 3


def test_cli_missing_file_exits_one(tmp_path, capsys):
    assert main(["plan", str(tmp_path / "nope.orbit")]) == 1
    assert "orbitkit:" in capsys.readouterr().err

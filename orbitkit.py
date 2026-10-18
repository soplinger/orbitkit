"""Offline two-body orbit toolkit.

Circular and escape speeds, Hohmann and bi-elliptic transfers, plane
changes, Kepler state solving, rocket-equation propellant budgets and
mission-file planning with line-numbered errors and CI exit codes.
Standard library only; no network access.
"""

import argparse
import json
import math
import sys

G0 = 9.80665  # standard gravity, m/s^2

# Gravitational parameter (km^3/s^2) and mean radius (km) per body.
BODIES = {
    "earth": {"mu": 398600.4418, "radius": 6378.137},
    "moon": {"mu": 4902.800066, "radius": 1737.4},
    "mars": {"mu": 42828.37362, "radius": 3396.2},
}

EXIT_OK = 0
EXIT_ERROR = 1
EXIT_BUDGET = 2


class MissionError(ValueError):
    """A single bad line in a mission file."""


def get_body(name):
    key = name.lower()
    if key not in BODIES:
        raise ValueError(f"unknown body '{name}' (expected one of {', '.join(sorted(BODIES))})")
    return BODIES[key]


def circular_speed(mu, r):
    """Speed (km/s) of a circular orbit of radius r (km)."""
    _positive(r, "orbit radius")
    return math.sqrt(mu / r)


def orbital_period(mu, a):
    """Period (s) of an orbit with semi-major axis a (km)."""
    _positive(a, "semi-major axis")
    return 2 * math.pi * math.sqrt(a ** 3 / mu)


def escape_speed(mu, r):
    """Escape speed (km/s) from radius r (km)."""
    _positive(r, "orbit radius")
    return math.sqrt(2 * mu / r)


def hohmann(mu, r1, r2):
    """Two-burn transfer between circular orbits r1 and r2 (km).

    Returns signed burns (km/s), their magnitude total and the transfer
    time (s). A negative burn means the burn points retrograde.
    """
    _positive(r1, "start radius")
    _positive(r2, "target radius")
    if math.isclose(r1, r2, rel_tol=1e-12):
        raise ValueError("start and target orbits are the same")
    a_t = (r1 + r2) / 2
    dv1 = math.sqrt(mu * (2 / r1 - 1 / a_t)) - math.sqrt(mu / r1)
    dv2 = math.sqrt(mu / r2) - math.sqrt(mu * (2 / r2 - 1 / a_t))
    return {
        "dv1_kms": dv1,
        "dv2_kms": dv2,
        "total_kms": abs(dv1) + abs(dv2),
        "time_s": orbital_period(mu, a_t) / 2,
    }


def bielliptic(mu, r1, r2, rb):
    """Three-burn transfer through an apex radius rb (km).

    Only worth using when rb is above both orbits, and the ratio of the
    target to the start radius is large (about 11.94 or more).
    """
    _positive(r1, "start radius")
    _positive(r2, "target radius")
    if rb <= max(r1, r2):
        raise ValueError("bi-elliptic apex must be above both orbits")
    a1 = (r1 + rb) / 2
    a2 = (r2 + rb) / 2
    dv1 = math.sqrt(mu * (2 / r1 - 1 / a1)) - math.sqrt(mu / r1)
    dv2 = math.sqrt(mu * (2 / rb - 1 / a2)) - math.sqrt(mu * (2 / rb - 1 / a1))
    dv3 = math.sqrt(mu / r2) - math.sqrt(mu * (2 / r2 - 1 / a2))
    return {
        "dv1_kms": dv1,
        "dv2_kms": dv2,
        "dv3_kms": dv3,
        "total_kms": abs(dv1) + abs(dv2) + abs(dv3),
        "time_s": math.pi * (math.sqrt(a1 ** 3 / mu) + math.sqrt(a2 ** 3 / mu)),
    }


def plane_change_dv(speed_kms, deg):
    """Delta-v (km/s) to rotate an orbit's plane by deg at a given speed."""
    if not 0 < deg <= 180:
        raise ValueError("plane change must be above 0 and at most 180 degrees")
    return 2 * speed_kms * math.sin(math.radians(deg) / 2)


def kepler_solve(mean_anomaly, e, tol=1e-12, max_iter=100):
    """Eccentric anomaly E (rad) from mean anomaly M (rad) by Newton's method."""
    if not 0 <= e < 1:
        raise ValueError("eccentricity must be at least 0 and below 1")
    m = math.fmod(mean_anomaly, 2 * math.pi)
    big_e = m if e < 0.8 else math.pi
    for _ in range(max_iter):
        step = (big_e - e * math.sin(big_e) - m) / (1 - e * math.cos(big_e))
        big_e -= step
        if abs(step) < tol:
            return big_e
    raise ValueError("Kepler solver did not converge")


def orbital_state(mu, a, e, mean_anomaly_deg):
    """Position and speed on an ellipse with semi-major axis a (km)."""
    _positive(a, "semi-major axis")
    big_e = kepler_solve(math.radians(mean_anomaly_deg), e)
    nu = 2 * math.atan2(
        math.sqrt(1 + e) * math.sin(big_e / 2),
        math.sqrt(1 - e) * math.cos(big_e / 2),
    )
    r = a * (1 - e * math.cos(big_e))
    v = math.sqrt(mu * (2 / r - 1 / a))
    return {
        "E_deg": math.degrees(big_e) % 360,
        "nu_deg": math.degrees(nu) % 360,
        "r_km": r,
        "v_kms": v,
    }


def rocket_propellant(dv_kms, isp_s, dry_kg):
    """Propellant mass (kg) for a burn of dv_kms using the rocket equation."""
    if dv_kms < 0:
        raise ValueError("delta-v cannot be negative")
    _positive(isp_s, "specific impulse")
    _positive(dry_kg, "dry mass")
    return dry_kg * (math.exp(dv_kms * 1000 / (G0 * isp_s)) - 1)


def plan_mission(text):
    """Parse and simulate a mission file.

    Returns a report dict. Errors are (line, message) pairs; when there
    are any, totals are left as None.
    """
    errors = []
    legs = []
    state = {
        "body": None,
        "body_name": None,
        "radius": None,
        "elapsed": 0.0,
        "isp": None,
        "dry": None,
        "budget": None,
    }

    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        parts = line.split()
        try:
            _apply_directive(parts, lineno, state, legs)
        except ValueError as exc:
            errors.append({"line": lineno, "message": str(exc)})

    if state["body"] is None and not errors:
        errors.append({"line": 0, "message": "mission has no 'body' directive"})

    report = {
        "body": state["body_name"],
        "legs": legs,
        "errors": errors,
        "total_kms": None,
        "elapsed_s": None,
        "propellant_kg": None,
        "budget_kms": state["budget"],
        "within_budget": None,
    }
    if errors:
        return report

    total = sum(leg["dv_kms"] for leg in legs)
    report["total_kms"] = total
    report["elapsed_s"] = state["elapsed"]
    if state["isp"] is not None:
        report["propellant_kg"] = rocket_propellant(total, state["isp"], state["dry"])
    if state["budget"] is not None:
        report["within_budget"] = total <= state["budget"] + 1e-9
    return report


def _apply_directive(parts, lineno, state, legs):
    verb = parts[0].lower()
    args = parts[1:]

    if verb == "body":
        _expect(args, 1, "body NAME")
        if legs or state["radius"] is not None:
            raise MissionError("'body' must come before any orbit or manoeuvre")
        state["body"] = get_body(args[0])
        state["body_name"] = args[0].lower()
        return

    if verb == "orbit":
        _expect(args, 2, "orbit circular ALT_KM")
        if args[0].lower() != "circular":
            raise MissionError("only 'orbit circular' is supported")
        if state["body"] is None:
            raise MissionError("set 'body' before 'orbit'")
        if legs:
            raise MissionError("'orbit' must come before any manoeuvre")
        alt = _number(args[1], "altitude")
        if alt < 0:
            raise MissionError("altitude cannot be negative")
        state["radius"] = state["body"]["radius"] + alt
        return

    if verb == "vehicle":
        _expect(args, 2, "vehicle ISP_S DRY_KG")
        state["isp"] = _positive(_number(args[0], "specific impulse"), "specific impulse")
        state["dry"] = _positive(_number(args[1], "dry mass"), "dry mass")
        return

    if verb == "budget":
        _expect(args, 1, "budget DV_KMS")
        state["budget"] = _non_negative(_number(args[0], "budget"), "budget")
        return

    if verb == "wait":
        _expect(args, 1, "wait HOURS")
        hours = _non_negative(_number(args[0], "hours"), "wait")
        state["elapsed"] += hours * 3600
        legs.append(_leg(lineno, f"wait {hours:g} h", 0.0, state))
        return

    if verb == "burn":
        _expect(args, 1, "burn DV_KMS")
        dv = _positive(_number(args[0], "delta-v"), "burn delta-v")
        legs.append(_leg(lineno, f"burn {dv:g} km/s", dv, state))
        return

    if verb == "hohmann":
        _require_orbit(state)
        _expect(args, 1, "hohmann TO_ALT_KM")
        target = _target_radius(args[0], state)
        try:
            result = hohmann(state["body"]["mu"], state["radius"], target)
        except ValueError as exc:
            raise MissionError(str(exc)) from None
        state["radius"] = target
        state["elapsed"] += result["time_s"]
        legs.append(_leg(lineno, f"hohmann to {args[0]} km", result["total_kms"], state))
        return

    if verb == "bielliptic":
        _require_orbit(state)
        _expect(args, 2, "bielliptic TO_ALT_KM VIA_ALT_KM")
        target = _target_radius(args[0], state)
        via = _target_radius(args[1], state)
        try:
            result = bielliptic(state["body"]["mu"], state["radius"], target, via)
        except ValueError as exc:
            raise MissionError(str(exc)) from None
        state["radius"] = target
        state["elapsed"] += result["time_s"]
        legs.append(_leg(lineno, f"bi-elliptic to {args[0]} km via {args[1]} km", result["total_kms"], state))
        return

    if verb == "planechange":
        _require_orbit(state)
        _expect(args, 1, "planechange DEG")
        deg = _number(args[0], "angle")
        speed = circular_speed(state["body"]["mu"], state["radius"])
        try:
            dv = plane_change_dv(speed, deg)
        except ValueError as exc:
            raise MissionError(str(exc)) from None
        legs.append(_leg(lineno, f"plane change {deg:g} deg", dv, state))
        return

    raise MissionError(f"unknown directive '{parts[0]}'")


def _leg(lineno, label, dv, state):
    return {"line": lineno, "label": label, "dv_kms": dv, "elapsed_s": state["elapsed"]}


def _require_orbit(state):
    if state["body"] is None:
        raise MissionError("set 'body' first")
    if state["radius"] is None:
        raise MissionError("set a starting 'orbit circular ALT_KM' first")


def _target_radius(token, state):
    alt = _non_negative(_number(token, "altitude"), "altitude")
    return state["body"]["radius"] + alt


def _expect(args, count, usage):
    if len(args) != count:
        raise MissionError(f"expected '{usage}'")


def _number(token, what):
    try:
        value = float(token)
    except ValueError:
        raise MissionError(f"{what} '{token}' is not a number") from None
    if not math.isfinite(value):
        raise MissionError(f"{what} must be finite")
    return value


def _positive(value, what):
    if value <= 0:
        raise ValueError(f"{what} must be greater than zero")
    return value


def _non_negative(value, what):
    if value < 0:
        raise ValueError(f"{what} cannot be negative")
    return value


def fmt_duration(seconds):
    hours = int(seconds // 3600)
    minutes = int(seconds % 3600 // 60)
    return f"{hours} h {minutes:02d} min"


def _cmd_info(args):
    body = get_body(args.body)
    r = body["radius"] + args.alt_km
    print(f"{args.body.title()} circular orbit at {args.alt_km:g} km")
    print(f"  speed:   {circular_speed(body['mu'], r):.3f} km/s")
    print(f"  period:  {orbital_period(body['mu'], r) / 60:.2f} min")
    print(f"  escape:  {escape_speed(body['mu'], r):.3f} km/s")
    return EXIT_OK


def _cmd_transfer(args):
    body = get_body(args.body)
    r1 = body["radius"] + args.from_alt
    r2 = body["radius"] + args.to_alt
    result = hohmann(body["mu"], r1, r2)
    print(f"Hohmann {args.from_alt:g} km -> {args.to_alt:g} km ({args.body.title()})")
    print(f"  burn 1:  {result['dv1_kms']:+.4f} km/s")
    print(f"  burn 2:  {result['dv2_kms']:+.4f} km/s")
    print(f"  total:   {result['total_kms']:.4f} km/s")
    print(f"  time:    {fmt_duration(result['time_s'])}")
    if args.via is not None:
        rb = body["radius"] + args.via
        bi = bielliptic(body["mu"], r1, r2, rb)
        print(f"Bi-elliptic via {args.via:g} km")
        print(f"  total:   {bi['total_kms']:.4f} km/s")
        print(f"  time:    {fmt_duration(bi['time_s'])}")
        winner = "bi-elliptic" if bi["total_kms"] < result["total_kms"] else "Hohmann"
        print(f"  cheaper: {winner}")
    return EXIT_OK


def _cmd_plane(args):
    dv = plane_change_dv(args.speed, args.deg)
    print(f"Plane change {args.deg:g} deg at {args.speed:g} km/s: {dv:.4f} km/s")
    return EXIT_OK


def _cmd_state(args):
    body = get_body(args.body)
    s = orbital_state(body["mu"], args.a_km, args.e, args.mean_anomaly)
    print(f"{args.body.title()} orbit a={args.a_km:g} km e={args.e:g} M={args.mean_anomaly:g} deg")
    print(f"  eccentric anomaly: {s['E_deg']:.4f} deg")
    print(f"  true anomaly:      {s['nu_deg']:.4f} deg")
    print(f"  radius:            {s['r_km']:.3f} km")
    print(f"  speed:             {s['v_kms']:.4f} km/s")
    return EXIT_OK


def _cmd_rocket(args):
    prop = rocket_propellant(args.dv, args.isp, args.dry)
    print(f"Burn {args.dv:g} km/s, Isp {args.isp:g} s, dry {args.dry:g} kg")
    print(f"  propellant: {prop:.1f} kg")
    print(f"  wet mass:   {args.dry + prop:.1f} kg")
    return EXIT_OK


def _cmd_plan(args):
    with open(args.file, encoding="utf-8") as handle:
        report = plan_mission(handle.read())

    if args.json:
        print(json.dumps(report, indent=2))
    else:
        _print_plan(report)

    if report["errors"]:
        return EXIT_ERROR
    if report["within_budget"] is False:
        return EXIT_BUDGET
    return EXIT_OK


def _print_plan(report):
    for err in report["errors"]:
        print(f"{err['line']}: error: {err['message']}")
    if report["errors"]:
        return
    print(f"Body: {report['body'].title()}")
    print(f"{'line':>5}  {'manoeuvre':<34} {'dv km/s':>9} {'total':>9} {'elapsed':>12}")
    running = 0.0
    for leg in report["legs"]:
        running += leg["dv_kms"]
        print(
            f"{leg['line']:>5}  {leg['label']:<34} {leg['dv_kms']:>9.4f} "
            f"{running:>9.4f} {fmt_duration(leg['elapsed_s']):>12}"
        )
    print(f"Total delta-v: {report['total_kms']:.4f} km/s")
    print(f"Elapsed:       {fmt_duration(report['elapsed_s'])}")
    if report["propellant_kg"] is not None:
        print(f"Propellant:    {report['propellant_kg']:.1f} kg")
    if report["budget_kms"] is not None:
        status = "within budget" if report["within_budget"] else "OVER BUDGET"
        print(f"Budget:        {report['budget_kms']:.4f} km/s ({status})")


def build_parser():
    parser = argparse.ArgumentParser(
        prog="orbitkit",
        description="Offline two-body orbit planner: speeds, transfers, plane changes, Kepler states and mission files.",
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("info", help="circular orbit speed, period and escape speed")
    p.add_argument("body", choices=sorted(BODIES))
    p.add_argument("alt_km", type=float, help="altitude above the surface, km")
    p.set_defaults(func=_cmd_info)

    p = sub.add_parser("transfer", help="Hohmann transfer, optionally vs bi-elliptic")
    p.add_argument("body", choices=sorted(BODIES))
    p.add_argument("from_alt", type=float, help="start altitude, km")
    p.add_argument("to_alt", type=float, help="target altitude, km")
    p.add_argument("--via", type=float, help="apex altitude for a bi-elliptic comparison, km")
    p.set_defaults(func=_cmd_transfer)

    p = sub.add_parser("plane", help="plane-change delta-v")
    p.add_argument("speed", type=float, help="orbital speed, km/s")
    p.add_argument("deg", type=float, help="plane change, degrees")
    p.set_defaults(func=_cmd_plane)

    p = sub.add_parser("state", help="position and speed from elements and mean anomaly")
    p.add_argument("body", choices=sorted(BODIES))
    p.add_argument("a_km", type=float, help="semi-major axis, km")
    p.add_argument("e", type=float, help="eccentricity, 0 to below 1")
    p.add_argument("mean_anomaly", type=float, help="mean anomaly, degrees")
    p.set_defaults(func=_cmd_state)

    p = sub.add_parser("rocket", help="propellant for a delta-v")
    p.add_argument("dv", type=float, help="delta-v, km/s")
    p.add_argument("isp", type=float, help="specific impulse, s")
    p.add_argument("dry", type=float, help="dry mass, kg")
    p.set_defaults(func=_cmd_rocket)

    p = sub.add_parser("plan", help="check and simulate a mission file")
    p.add_argument("file", help="mission file")
    p.add_argument("--json", action="store_true", help="print the report as JSON")
    p.set_defaults(func=_cmd_plan)

    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except (ValueError, OSError) as exc:
        print(f"orbitkit: {exc}", file=sys.stderr)
        return EXIT_ERROR


if __name__ == "__main__":
    sys.exit(main())

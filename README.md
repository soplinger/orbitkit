# Orbitkit

**Offline two-body orbit planner: speeds, transfers, plane changes, Kepler states and mission files, all from one stdlib Python file.**

Orbitkit answers the everyday questions of orbital mechanics without a network connection or a physics package. Ask it for a circular orbit's speed, the cost of a Hohmann or bi-elliptic transfer, or the propellant a burn needs. Or write a short plain-text mission file and get the full delta-v budget, elapsed time and line-numbered errors back. It uses only the Python standard library.

## Features

- **Circular orbit figures**: speed, period and escape speed at any altitude for Earth, the Moon and Mars.
- **Hohmann transfers**: two-burn delta-v, signed burns (prograde or retrograde) and transfer time between any two circular orbits.
- **Bi-elliptic comparison**: three-burn transfer through an apex radius, with an automatic "which is cheaper" verdict.
- **Plane-change delta-v** at a given orbital speed, with angle validation.
- **Kepler state solver**: eccentric anomaly by Newton's method, then true anomaly, radius and speed for any ellipse.
- **Rocket equation**: propellant and wet mass for a delta-v, given specific impulse and dry mass.
- **Mission files**: plain-text manoeuvre lists (`body`, `orbit`, `hohmann`, `bielliptic`, `planechange`, `burn`, `wait`, `vehicle`, `budget`).
- **Simulated timeline**: each manoeuvre's delta-v, running total and elapsed time, printed as a table.
- **Budget gate**: flags a mission whose total delta-v exceeds the `budget` directive.
- **Line-numbered lint**: every bad line is reported with `line: error: message`, and totals are withheld until the file is clean.
- **JSON output** for `plan --json`, for scripts and dashboards.
- **CI exit codes**: `0` clean, `1` file or lint error, `2` over the delta-v budget.
- **Comments everywhere**: `#` starts a comment, whole-line or inline.

## Quick start

```sh
python3 orbitkit.py info earth 200
python3 orbitkit.py transfer earth 200 35786 --via 100000
python3 orbitkit.py plan examples/leo_to_geo.orbit
```

The example mission (`examples/leo_to_geo.orbit`) is a 200 km parking orbit, a Hohmann transfer to geostationary altitude, a two-hour coast, a 0.1° plane trim, and a propellant estimate for a 320 s engine with a 1500 kg dry mass, against a 4.5 km/s budget:

```text
body earth
orbit circular 200        # 200 km parking orbit
hohmann 35786             # to geostationary altitude
wait 2                    # coast before the trim
planechange 0.1           # small plane correction
vehicle 320 1500          # Isp in seconds, dry mass in kg
budget 4.5                # mission delta-v limit, km/s
```

`plan` prints one row per manoeuvre, then the totals, the propellant and the budget status. The Hohmann leg is about 3.93 km/s and 5 h 15 min, and the whole mission is about 3.94 km/s over 7 h 15 min. Run the command to see the exact figures.

## Mission file reference

| Directive | Meaning |
| --- | --- |
| `body NAME` | `earth`, `moon` or `mars`. Must come first. |
| `orbit circular ALT_KM` | Starting parking orbit. Must come before any manoeuvre. |
| `hohmann TO_ALT_KM` | Two-burn transfer to a new circular orbit. |
| `bielliptic TO_ALT_KM VIA_ALT_KM` | Three-burn transfer through an apex above both orbits. |
| `planechange DEG` | Plane change at the current orbit's speed. |
| `burn DV_KMS` | Any other impulsive burn, counted at its magnitude. |
| `wait HOURS` | Coast time, which adds to elapsed time. |
| `vehicle ISP_S DRY_KG` | Enables propellant estimates. |
| `budget DV_KMS` | Enables the CI budget gate. |

## Commands

| Command | Does |
| --- | --- |
| `info BODY ALT_KM` | Circular speed, period and escape speed. |
| `transfer BODY FROM_ALT TO_ALT [--via ALT]` | Hohmann transfer, plus bi-elliptic comparison with `--via`. |
| `plane SPEED_KMS DEG` | Plane-change delta-v. |
| `state BODY A_KM E M_DEG` | Position and speed from semi-major axis, eccentricity and mean anomaly. |
| `rocket DV_KMS ISP DRY_KG` | Propellant and wet mass. |
| `plan FILE [--json]` | Lint and simulate a mission file. |

## How to run

Requires Python 3.8 or later. There are no dependencies to install.

```sh
python3 orbitkit.py --help
python3 orbitkit.py state earth 26561.75 0.7 30
python3 orbitkit.py rocket 3.9 320 1500
python3 orbitkit.py plan examples/leo_to_geo.orbit --json
```

To gate a CI job on a mission file, run `plan` and check the exit code. Exit `1` means the file has errors, and `2` means the total delta-v is over the `budget` line.

## How to test

Tests use pytest and live in `tests/`:

```sh
python3 -m pytest
```

The suite covers the orbital formulas, Kepler solving, the rocket equation, mission parsing and line-numbered errors, the budget gate, and the command-line exit codes.

## Project layout

```text
orbitkit.py              the whole tool: pure functions plus the argparse CLI
examples/leo_to_geo.orbit  a worked mission file
tests/                   pytest suite
```

## Notes and limits

- Two-body, impulsive burns only: no atmospheric drag, no perturbations, no finite-burn losses.
- Altitudes are above the body's mean radius, and all distances are in kilometres.
- The Hohmann and bi-elliptic transfers are coplanar. Put plane changes in their own `planechange` line.

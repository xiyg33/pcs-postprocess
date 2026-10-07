"""Public SIM/HIL comparison command; each mode owns its metrics."""

import argparse
from datetime import datetime
import json
from pathlib import Path


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--mode", default="pwr", choices=("pwr", "frt", "freq_reg", "inertia", "volt_prot"))
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--case", type=int)
    selection.add_argument("--all", action="store_true")
    parser.add_argument("--paths-config", type=Path)
    parser.add_argument("--output", type=Path)


def execute(args, parser: argparse.ArgumentParser) -> int:
    if args.mode == "frt":
        from . import compare_frt as mode
    elif args.mode == "freq_reg":
        from . import compare_freq_reg as mode
    elif args.mode == "inertia":
        from . import compare_inertia as mode
    elif args.mode == "volt_prot":
        from . import compare_volt_prot as mode
    else:
        from . import compare_pwr as mode
    output = args.output or Path.cwd() / "result_compare" / args.mode / datetime.now().strftime("%Y-%m-%d_%H%M%S_%f")
    try:
        result = mode.run(args.paths_config, output, args.case)
    except (ValueError, KeyError, OSError) as exc:
        parser.exit(2, f"error: {exc}\n")
    print(json.dumps({"output": str(output.resolve()), "completed": result["completed"],
                      "skipped": result.get("skipped", 0), "exceptions": result["exceptions"]}))
    return int(bool(result["exceptions"]))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    return execute(parser.parse_args(argv), parser)


if __name__ == "__main__":
    raise SystemExit(main())

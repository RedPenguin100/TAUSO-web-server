"""Check that this web server and the TAUSO it is pinned to still fit together.

Run it inside the image, after a build and before trusting a deploy:

    docker compose exec -T tauso-web /opt/conda/bin/python /app/smoke_test.py

It designs a handful of ASOs against a small target, which is enough to catch the
things that go wrong between the two: an API that was renamed, a model asking for
features the pipeline no longer produces, a data file whose schema moved, a native
binary missing from the image. Each of those is silent until a real job hits it,
and a real job costs minutes and an email to a user.
"""

import shutil
import sys
import traceback

CHECKS = []


def check(name):
    def register(fn):
        CHECKS.append((name, fn))
        return fn

    return register


@check("tauso API the pipeline imports")
def _api():
    from tauso.aso_generation import default_config, design_asos, summarize_design, tox_details

    for fn in (default_config, design_asos, summarize_design, tox_details):
        assert callable(fn), fn
    return "design_asos, summarize_design, tox_details, default_config"


@check("bowtie on PATH")
def _bowtie():
    found = shutil.which("bowtie")
    assert found, "bowtie is not on PATH; the off-target search will fail at run time"
    return found


@check("model present and its features computable")
def _model():
    from tauso.inference import load_model

    booster, features = load_model()
    assert features, "the model reports no features"
    return f"{len(features)} features"


@check("a design runs end to end")
def _design():
    from tauso.aso_generation import design_asos, summarize_design, tox_details

    import pipeline_runner as runner

    config = runner.JobConfig(
        target_data="", target_mrna_name="MALAT1", source_info="smoke test", user_email="smoke@local"
    )
    design_config = _design_config(runner, config)
    ranked = design_asos(
        config.target_mrna_name,
        cell_line="A549",
        aso_sizes=[len(config.chemical_pattern)],
        config=design_config,
        first_n=5,
        off_targets=False,
    )
    assert len(ranked) == 5, f"expected 5 candidates, got {len(ranked)}"
    designed, safety = summarize_design(ranked), tox_details(ranked)
    assert len(designed) == 5 and len(safety) == 5
    score = designed[designed.columns[-1]]
    assert score.notna().all(), "some candidates scored NaN"
    return f"{len(designed)} scored, {len(designed.columns)} + {len(safety.columns)} result columns"


@check("blank conditions reach the model as missing")
def _blank_conditions():
    from tauso.aso_generation import design_asos

    import pipeline_runner as runner

    config = runner.JobConfig(
        target_data="", target_mrna_name="MALAT1", source_info="smoke test", user_email="smoke@local"
    )
    ranked = design_asos(
        config.target_mrna_name,
        cell_line="None",
        aso_sizes=[len(config.chemical_pattern)],
        config=_design_config(runner, config),
        first_n=5,
        off_targets=False,
    )
    missing = ["volume_nm", "density_cells_per_well", "transfection_gymnosis"]
    for column in missing:
        assert ranked[column].isna().all(), f"{column} should be missing when unset"
    return ", ".join(missing)


@check("sugar and backbone validation")
def _patterns():
    from pipeline_runner import CHEMISTRIES, describe_pattern_problem

    for name, spec in CHEMISTRIES.items():
        problem = describe_pattern_problem(spec["pattern"], spec["ps_pattern"])
        assert problem is None, f"{name} is rejected: {problem}"
    assert describe_pattern_problem("MMMMMddddddddddMMMMM", "*" * 20), "a wrong-length backbone passed"
    assert describe_pattern_problem("MMMMMMMMMMMMMMMMMMMM", "*" * 19), "a non-gapmer passed"
    return f"{len(CHEMISTRIES)} chemistries accepted, bad patterns rejected"


def _design_config(runner, config):
    design_config = runner.default_config()
    design_config.standard_chemical_pattern = config.chemical_pattern
    design_config.standard_ps_pattern = config.ps_pattern
    design_config.standard_modification = config.modification
    design_config.transfection_method = config.transfection
    design_config.volume = float("nan") if config.dosage_nm is None else config.dosage_nm
    design_config.cell_per_well = float("nan") if config.cell_density is None else config.cell_density
    return design_config


def main():
    sys.path.insert(0, "/app")
    failures = 0
    for name, fn in CHECKS:
        try:
            print(f"  {name} ... {fn()}", flush=True)
        except Exception as exc:
            failures += 1
            print(f"  {name} ... FAILED: {type(exc).__name__}: {exc}", flush=True)
            traceback.print_exc()
    print(f"\n{len(CHECKS) - failures}/{len(CHECKS)} checks passed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())

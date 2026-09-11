"""Missing directions and verifier unknowns must not become benchmark failures."""
import contextlib
import io
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bench'))
import golden_tolerance_benchmark as benchmark


class GoldenToleranceAccountingTests(unittest.TestCase):
    def test_partial_searches_and_unknown_queries_keep_separate_denominators(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run, archived = root / 'run', root / 'archived'

            def result(variant, index, direction, recovery, cpu):
                folder = benchmark.result_folder(variant, dict(dataset='golden', index=index, direction=direction))
                benchmark.save(folder / 'search.json', dict(complete=True, rows=[dict(
                    compute_cpu_excluding_persistence_and_loading_seconds=cpu, capped=False)]))
                benchmark.save(folder / 'full_sweep_evaluation.json', dict(reference_recovery=recovery))

            control, trial = [run / 'runs' / k for k in benchmark.TOLERANCES]
            # A recovered direction proves union recovery despite a missing
            # reverse search, but this case cannot enter paired CPU totals.
            result(control, 0, 'R_to_P', 'recovered', 50)
            for direction in ('R_to_P', 'P_to_R'):
                result(trial, 0, direction, 'not_recovered', 100)
            # A completed search with an unknown verifier result remains
            # eligible for CPU comparison, while coverage stays unknown.
            for direction in ('R_to_P', 'P_to_R'):
                result(control, 1, direction, 'not_recovered', 2)
                result(trial, 1, direction, 'unknown' if direction == 'R_to_P' else 'not_recovered', 3)
            with patch.object(benchmark, 'BASELINE', archived), contextlib.redirect_stdout(io.StringIO()):
                benchmark.collect(SimpleNamespace(run=run))
            summary = benchmark.read(run / 'comparison/summary.json')
            rows = benchmark.read(run / 'comparison/case_metrics.json')
            self.assertEqual(summary['paired_complete_cases'], 1)
            self.assertEqual(summary['paired_mapping_cpu'], {'tol_1p0': 4, 'tol_1p5': 6})
            self.assertEqual(rows[0]['variants']['tol_1p0']['recovery'], 'recovered')
            self.assertIsNone(rows[0]['variants']['tol_1p0']['mapping_cpu'])
            self.assertEqual(rows[1]['variants']['tol_1p5']['recovery'], 'unknown')
            self.assertEqual(rows[2]['variants']['tol_1p5']['recovery'], 'unknown')
            self.assertEqual(summary['metrics']['tol_1p5']['recovery']['not_recovered'], 1)


if __name__ == '__main__':
    unittest.main()

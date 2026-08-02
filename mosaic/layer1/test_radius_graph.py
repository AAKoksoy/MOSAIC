"""Regression tests for Layer 1 spatial graph summaries."""

import pandas as pd
import unittest

from mosaic.layer1.radius_graph import summarize_spatial_graph


class TestSpatialGraphSummary(unittest.TestCase):
    """Verify summary statistics for sample-specific graphs."""

    def test_degree_counts_are_calculated_per_sample(self) -> None:
        """Source and target indices align without mixing repeated IDs."""

        cells = pd.DataFrame(
            {
                "sample_id": ["sample_a"] * 4 + ["sample_b"] * 3,
                "cell_id": ["cell_1", "cell_2", "cell_3", "cell_4"]
                + ["cell_1", "cell_2", "cell_3"],
            }
        )
        edges = pd.DataFrame(
            {
                "sample_id": ["sample_a", "sample_a", "sample_b"],
                "source": ["cell_1", "cell_2", "cell_1"],
                "target": ["cell_2", "cell_3", "cell_3"],
                "distance_um": [10.0, 12.0, 8.0],
                "relationship_type": ["NEAR_TO"] * 3,
            }
        )

        summary = summarize_spatial_graph(cells, edges)

        # Per-sample degrees are [1, 2, 1, 0] and [1, 0, 1].
        self.assertEqual(summary["number_of_nodes"], 7)
        self.assertEqual(summary["number_of_edges"], 3)
        self.assertEqual(summary["number_of_samples"], 2)
        self.assertEqual(summary["number_of_isolated_nodes"], 2)
        self.assertAlmostEqual(summary["fraction_isolated"], 2 / 7)
        self.assertEqual(summary["minimum_degree"], 0)
        self.assertEqual(summary["maximum_degree"], 2)
        self.assertAlmostEqual(summary["mean_degree"], 6 / 7)
        self.assertEqual(summary["median_degree"], 1)


if __name__ == "__main__":
    unittest.main()

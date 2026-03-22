import os
import tempfile
import tomllib
import unittest
from pathlib import Path

from aloha.scenario import Scenario

SCENARIOS_DIR = Path(__file__).parent.parent.parent / "scenarios"
TOML_FILES = list(SCENARIOS_DIR.glob("*.toml"))


class TestScenario(unittest.TestCase):
    def test_scenario_constructor_no_args(self):
        scenario = Scenario()
        self.assertEqual(scenario.scenario, {})

    def test_scenario_constructor_dict(self):
        data = {"options": {"test": True}}
        scenario = Scenario(data)
        self.assertEqual(scenario.scenario, data)

    def test_scenario_constructor_str(self):
        for scenario_file in TOML_FILES:
            with self.subTest(file=scenario_file.name):
                scenario = Scenario(str(scenario_file))
                self.assertTrue(scenario.scenario)

    def test_scenario_constructor_path(self):
        for scenario_file in TOML_FILES:
            with self.subTest(file=scenario_file.name):
                scenario = Scenario(scenario_file)
                self.assertTrue(scenario.scenario)

    def test_scenario_constructor_invalid(self):
        with self.assertRaises(ValueError):
            Scenario(123)


if __name__ == "__main__":
    unittest.main()

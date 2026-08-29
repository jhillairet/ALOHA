import os
import tempfile
import unittest
from pathlib import Path

# Use tomllib for Python 3.11+, tomli for earlier versions
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

from aloha.scenario import Scenario

SCENARIOS_DIR = Path(__file__).parent.parent.parent / "scenarios"
TOML_SCENARIO_FILES = list(SCENARIOS_DIR.glob("*.toml"))


class TestScenario(unittest.TestCase):
    def test_scenario_constructor_no_args(self):
        scenario = Scenario()
        self.assertEqual(scenario.scenario, {})

    def test_scenario_constructor_dict(self):
        data = {"options": {"test": True}}
        scenario = Scenario(data)
        self.assertEqual(scenario.scenario, data)

    def test_scenario_constructor_str(self):
        for scenario_file in TOML_SCENARIO_FILES:
            with self.subTest(file=scenario_file.name):
                scenario = Scenario(str(scenario_file))
                self.assertTrue(scenario.scenario)

    def test_scenario_constructor_path(self):
        for scenario_file in TOML_SCENARIO_FILES:
            with self.subTest(file=scenario_file.name):
                scenario = Scenario(scenario_file)
                self.assertTrue(scenario.scenario)

    def test_scenario_constructor_invalid(self):
        with self.assertRaises(ValueError):
            Scenario(123)


if __name__ == "__main__":
    unittest.main()

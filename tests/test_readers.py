# -*- coding: utf-8 -*-
"""
unit tests for reading different file formats in DATA/<program>/
"""
import unittest
import logging

# # Logging
logger = logging.getLogger(__name__)
logging.basicConfig(format="[testing] %(message)s", level=logging.ERROR)

from sgdml_dataset_generation.readers.qchem import QChemOutputFile

class TestQChemOutputFile(unittest.TestCase):
    def test_parser(self):
        log_files = ["DATA/QChem/methylium.out", "DATA/QChem/squaraine.out"]
        for log_file in log_files:
            logger.info(f"reading {log_file}")
            with open(log_file) as f:
                results = QChemOutputFile(f)
            # check that required keys exist
            results['geometry (au)']
            results['energies (au)']
            results['gradients (au)']
            results['derivative couplings (au)']
                
if __name__ == "__main__":
    unittest.main()

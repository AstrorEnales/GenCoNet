#!/usr/bin/env python3

import unittest
from utils import name_utils


class TestMethods(unittest.TestCase):
    def test_normalize_drug_names(self):
        self.assertEqual(name_utils.normalize_drug_names({'Dicumarol'}), {'Dicumarol'})
        self.assertEqual(name_utils.normalize_drug_names({'Dicumarol', 'dicumarol'}), {'Dicumarol'})
        self.assertEqual(name_utils.normalize_drug_names({'Phenolphthalein', 'PHENOLPHTHALEIN'}), {'Phenolphthalein'})
        self.assertEqual(name_utils.normalize_drug_names(
            {'fludrocortisone acetate', 'Fludrocortisone acetate', 'Fludrocortisone Acetate'}),
            {'Fludrocortisone Acetate'})

    def test_drug_names_synonym(self):
        self.assertFalse(name_utils.drug_names_synonym({'Dicumarol', 'Xanthinol'}))
        self.assertTrue(name_utils.drug_names_synonym({'Xanthinol', 'xanthinol'}))
        self.assertTrue(name_utils.drug_names_synonym({'Xantinol', 'xanthinol'}))
        self.assertTrue(name_utils.drug_names_synonym({'Carphenazine', 'carfenazine'}))
        self.assertTrue(name_utils.drug_names_synonym({'aciclovir', 'Acyclovir'}))
        self.assertTrue(name_utils.drug_names_synonym({'metacycline', 'Methacycline'}))
        self.assertTrue(name_utils.drug_names_synonym({'debrisoquine', 'Debrisoquin'}))
        self.assertTrue(name_utils.drug_names_synonym({'Aspergillus niger allergenic extract',
                                                       'ALLERGENIC EXTRACT, ASPERGILLUS NIGER'}))
        self.assertTrue(name_utils.drug_names_synonym({'Sodium Phosphate, Monobasic',
                                                       'SODIUM PHOSPHATE,MONOBASIC'}))


if __name__ == '__main__':
    unittest.main()

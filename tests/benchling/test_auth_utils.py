import unittest
from unittest.mock import patch


class TestArchiveOligo(unittest.TestCase):
    def setUp(self):
        self.example = "example"

    def test_example(self):
        # Arrange
        starting_data = self.example
        # Act
        acted_data = example(starting_data)
        # Assert
        self.assertDictEqual(acted_data, self.example_acted_data)


if __name__ == '__main__':
    unittest.main()
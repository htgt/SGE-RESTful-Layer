import unittest
from unittest.mock import patch
from src.benchling.archive_entity import (
    archive_entity, 
    archive_oligo, 
    send_archive_request, 
    prepare_archive_json, 
    generate_entity_ids_str
)


class TestArchiveOligo(unittest.TestCase):
    def setUp(self):
        self.example = "example"

    def test_archive_oligo(self):
        # Arrange
        starting_data = self.example
        # Act
        acted_data = example(starting_data)
        # Assert
        self.assertDictEqual(acted_data, self.example_acted_data)
        
    def test_archive_entity(self):
        # Arrange
        starting_data = self.example
        # Act
        acted_data = example(starting_data)
        # Assert
        self.assertDictEqual(acted_data, self.example_acted_data)
    
    def test_send_archive_request(self):
        # Arrange
        starting_data = self.example
        # Act
        acted_data = example(starting_data)
        # Assert
        self.assertDictEqual(acted_data, self.example_acted_data)
    
    def test_prepare_archive_json(self):
        # Arrange
        starting_data = self.example
        # Act
        acted_data = example(starting_data)
        # Assert
        self.assertDictEqual(acted_data, self.example_acted_data)
            
    def test_generate_entity_ids_str(self):
        # Arrange
        starting_data = self.example
        # Act
        acted_data = example(starting_data)
        # Assert
        self.assertDictEqual(acted_data, self.example_acted_data)


if __name__ == '__main__':
    unittest.main()

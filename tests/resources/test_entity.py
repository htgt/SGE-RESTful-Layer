import unittest

from mock import patch

from src.resources.entity import Entity


class TestEntity(unittest.TestCase):
    def setUp(self):
        return

    @patch('src.resources.entity.entities', [{'id': 1}])
    def test_get_200(self):
        # arrange
        test_input200 = 1
        expected = ({'entity': {'id': 1}}, 200)

        entity = Entity()

        # act
        actual200 = entity.get(test_input200)

        # assert
        self.assertEqual(actual200, expected)

    def test_get_404(self):
        # arrange
        test_input404 = 888
        expected = ({'entity': None}, 404)

        entity = Entity()

        # act
        actual404 = entity.get(test_input404)

        # assert
        self.assertEqual(actual404, expected)


if __name__ == '__main__':
    unittest.main()
